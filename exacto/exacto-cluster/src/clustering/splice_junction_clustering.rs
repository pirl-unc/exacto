// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use bimap::BiMap;
use exacto_core::prelude::{overlaps, FastaMap, Strand};
use exacto_caller::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::clustering::junction_placement::{place_read_junctions, pool_insertion_abutting_junctions, PlacedJunction};
use crate::clustering::splice_transition::SpliceTransitionReadCounts;
use crate::reference::junction_orientation::{is_aligned_against_transcript, mirror_splice_junction};
use crate::reference::splice_junction_annotation_index::SpliceJunctionAnnotationIndex;
use crate::clustering::splice_junction_cluster::SpliceJunctionCluster;


pub fn cluster_transcript_models_by_splice_junctions(
    transcript_models: &Vec<TranscriptModel>,
    annotation_index: &SpliceJunctionAnnotationIndex,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    min_reads_per_cluster: usize,
    expected_sequencing_error: f64,
    max_fpr: f64,
    max_unspliced_locus_gap: u32,
    dna_carrier_read_ids: Option<&HashSet<usize>>
) -> Vec<SpliceJunctionCluster> {
    // Step 1. Build chains: each read's junctions in transcript order.
    // HashMap<junction chain, HashSet<read ID>>
    let mut chains_to_reads: HashMap<Vec<SpliceJunction>, HashSet<usize>> = HashMap::new();
    let mut cycle_junctions: HashSet<SpliceJunction> = HashSet::new();
    let mut chained_read_ids: HashSet<usize> = HashSet::new();
    let mut placed_junctions: Vec<Vec<PlacedJunction>> = transcript_models
        .iter()
        .map(|transcript_model| place_read_junctions(transcript_model, chromosome_names_map, fasta_map))
        .collect();
    pool_insertion_abutting_junctions(&mut placed_junctions);

    // The reference genes of every annotated junction and the annotated junction pairs.
    let junctions_to_genes: &HashMap<(u16, u32, u32), BTreeSet<Box<str>>> = &annotation_index.junctions_to_genes;
    let annotated_transitions: &HashSet<(SpliceJunction, SpliceJunction)> = &annotation_index.annotated_transitions;

    // Hold every novel splice junction to its minimum read support, on the placements the
    // pooling settled. A read splicing one that falls short is removed: it keys no chain and
    // it is not a junction-less read.
    let read_splice_junctions: Vec<Vec<SpliceJunction>> = placed_junctions
        .iter()
        .map(|junctions| junctions.iter().map(|placed| placed.junction.clone()).collect())
        .collect();
    let novel_junction_read_counts: RNANovelJunctionReadCounts = RNANovelJunctionReadCounts::new(
        &read_splice_junctions,
        |intron| junctions_to_genes.contains_key(intron)
    );
    let novel_junction_read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &novel_junction_read_counts.get_depths(),
        expected_sequencing_error,
        max_fpr,
        1
    );
    let novel_junction_filter: RNAJunctionReadSupportFilter<RNANovelJunctionReadCounts> = RNAJunctionReadSupportFilter::new(
        &novel_junction_read_counts,
        &novel_junction_read_support_index
    );
    let unsupported_reads: HashSet<usize> = read_splice_junctions
        .iter()
        .enumerate()
        .filter(|(_, splice_junctions)| !novel_junction_filter.passes(splice_junctions))
        .map(|(index, _)| index)
        .collect();

    // The back-splices of every read, at their read positions.
    // Vec<Vec<(read position, back-splice)>>
    let back_splices: Vec<Vec<(u32, SpliceJunction)>> = transcript_models
        .iter()
        .map(|transcript_model| {
            transcript_model
                .get_variant_records()
                .iter()
                .filter(|variant_record| {
                    *variant_record.get_variant_type() == VariantType::CircularRNA
                        && variant_record.get_chromosome_1() == variant_record.get_chromosome_2()
                })
                .map(|variant_record| {
                    let orient = |strand: &Strand| -> Strand {
                        match (variant_record.get_strand_1(), strand) {
                            (Strand::Reverse, Strand::Forward) => Strand::Reverse,
                            (Strand::Reverse, Strand::Reverse) => Strand::Forward,
                            _ => strand.clone()
                        }
                    };
                    let cycle: SpliceJunction = SpliceJunction::new(
                        variant_record.get_chromosome_1(),
                        variant_record.get_chromosome_2(),
                        variant_record.get_position_1(),
                        variant_record.get_position_2(),
                        orient(variant_record.get_strand_1()),
                        orient(variant_record.get_strand_2())
                    );
                    (variant_record.get_read_position_1(), cycle)
                })
                .collect()
        })
        .collect();

    for (index, ((transcript_model, junctions), cycles)) in transcript_models
        .iter()
        .zip(placed_junctions.into_iter())
        .zip(back_splices.into_iter())
        .enumerate()
    {
        if unsupported_reads.contains(&index) {
            chained_read_ids.insert(transcript_model.get_read_id());
            continue;
        }

        // (read position, junction number, junction). The transcript model numbers its
        // junctions by read position, so read position is transcript order, and it is the
        // one order a back-splice can be slotted into.
        let mut positioned: Vec<(u32, Option<u16>, SpliceJunction)> = junctions
            .into_iter()
            .map(|placed| (placed.read_position, Some(placed.number), placed.junction))
            .collect();
        for (read_position, cycle) in cycles.into_iter() {
            positioned.push((read_position, None, cycle));
        }
        if positioned.is_empty() {
            continue;   // junction-less reads are keyed in Step 9
        }
        positioned.sort_by_key(|(read_position, _, _)| *read_position);
        debug_assert!(
            positioned
                .iter()
                .filter_map(|(_, number, _)| *number)
                .enumerate()
                .all(|(i, number)| number as usize == i + 1),
            "Splice junction numbers of read {} are not 1..=n in read order.",
            transcript_model.get_read_id()
        );
        // Key the chain in transcript orientation, so that the reads of one transcript sequenced
        // from either end share it: a read aligned against the strand its splice motifs give
        // has its chain reversed and every intron entered from the other end. A back-splice key
        // is already the same from both ends.
        let introns: Vec<SpliceJunction> = positioned
            .iter()
            .filter(|(_, number, _)| number.is_some())
            .map(|(_, _, junction)| junction.clone())
            .collect();
        let mut numbered: Vec<(Option<u16>, SpliceJunction)> = positioned
            .into_iter()
            .map(|(_, number, junction)| (number, junction))
            .collect();
        if is_aligned_against_transcript(&introns, chromosome_names_map, fasta_map) {
            numbered = numbered
                .iter()
                .rev()
                .map(|(number, junction)| match number {
                    Some(_) => (*number, mirror_splice_junction(junction)),
                    None => (None, junction.clone())
                })
                .collect();
        }
        let chain: Vec<SpliceJunction> = numbered
            .into_iter()
            .map(|(number, junction)| {
                if number.is_none() {
                    cycle_junctions.insert(junction.clone());
                }
                junction
            })
            .collect();
        chained_read_ids.insert(transcript_model.get_read_id());
        chains_to_reads
            .entry(chain)
            .or_default()
            .insert(transcript_model.get_read_id());
    }

    // Hold every ordered pair of consecutive junctions that a chain spells to its minimum read
    // support. A chain spelling one that falls short is dropped.
    let transition_read_counts: SpliceTransitionReadCounts = SpliceTransitionReadCounts::new(
        &chains_to_reads,
        annotated_transitions,
        junctions_to_genes
    );
    let transition_read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &transition_read_counts.get_depths(),
        expected_sequencing_error,
        max_fpr,
        1
    );
    let transition_filter: RNAJunctionReadSupportFilter<SpliceTransitionReadCounts> = RNAJunctionReadSupportFilter::new(
        &transition_read_counts,
        &transition_read_support_index
    );
    chains_to_reads.retain(|chain, _| transition_filter.passes(chain));

    let chains: HashSet<&Vec<SpliceJunction>> = chains_to_reads.keys().collect();

    // Step 2. Build a map of junctions to chains.
    // HashMap<junction, HashSet<junction chain>>
    let mut junctions_to_chains: HashMap<&SpliceJunction, HashSet<&Vec<SpliceJunction>>> = HashMap::new();
    for &chain in chains.iter() {
        for junction in chain.iter() {
            junctions_to_chains
                .entry(junction)
                .or_default()
                .insert(chain);
        }
    }

    // Step 3. Identify superchains for each chain.
    //
    // A superchain of chain C is any strictly-longer chain that contains C as a *contiguous*
    // subchain. It is a fuller isoform that C could be a 5'/3'-truncated fragment of. A chain with
    // no superchains of its own is *maximal* (identified in Step 4): a full-length backbone
    // contained in nothing longer. Example (letters are junctions in transcript order):
    //
    //   Chain   Junctions       Superchains of it   Maximal?
    //   C1      [B, C]          C2, C3, C4          no
    //   C2      [A, B, C]       C3                  no
    //   C4      [B, C, D]       C3                  no
    //   C3      [A, B, C, D]    — (none)            yes
    //
    // C2 and C4 are superchains of C1 but are not maximal, because C3 contains them; only C3
    // is maximal. Containment is transitive, so every non-maximal chain has at least one
    // maximal superchain — the guarantee Step 6 relies on.

    /// True if `a` is a strictly shorter, contiguous run inside `b`.
    /// This is the degradation model: a 5' or 3' truncated read keeps a contiguous block of
    /// its parent's junctions (a prefix, suffix, or internal window), so a degradation fragment
    /// is a contiguous subchain of the full-length isoform.
    /// A non-contiguous subset (e.g. a missing middle junction) stays its own chain.
    fn is_contiguous_subchain(a: &[SpliceJunction], b: &[SpliceJunction]) -> bool {
        if a.len() >= b.len() {
            return false;
        }
        let a_length: usize = a.len();
        let b_length: usize = b.len();
        for i in 0..(b_length - a_length + 1) {
            if b[i..i + a_length] == a[..] {
                return true;
            }
        }
        false
    }

    /// Return every chain that contains `chain` as a *contiguous* subchain (i.e. superchains).
    fn identify_superchains<'a>(
        chain: &'a Vec<SpliceJunction>,
        junctions_to_chains: &HashMap<&'a SpliceJunction, HashSet<&'a Vec<SpliceJunction>>>
    ) -> HashSet<&'a Vec<SpliceJunction>> {
        let mut containers: HashSet<&'a Vec<SpliceJunction>> = HashSet::new();
        for signatures in junctions_to_chains.get(chain.first().unwrap()) {
            for &signature_ in signatures {
                if is_contiguous_subchain(chain, signature_) {
                    containers.insert(signature_);
                }
            }
        }
        containers
    }

    // HashMap<junction chain, HashSet<superchain>>
    let mut chains_to_superchains: HashMap<&Vec<SpliceJunction>, HashSet<&Vec<SpliceJunction>>> = HashMap::new();
    for &chain in chains.iter() {
        chains_to_superchains.insert(chain, identify_superchains(chain, &junctions_to_chains));
    }

    // Step 4. Identify maximal superchains (full-length isoform backbones - truncations of nothing).
    // These are the chains with no superchains of their own.
    let mut maximal_superchains: HashSet<&Vec<SpliceJunction>> = HashSet::new();
    for &chain in chains.iter() {
        if chains_to_superchains[chain].is_empty() {
            maximal_superchains.insert(chain);
        }
    }

    // Step 5. Identify the reference genes of each chain, from the map of junctions to
    // reference genes.

    /// The annotated genes a chain's junctions belong to.
    fn chain_genes<'a>(
        chain: &[SpliceJunction],
        junctions_to_genes: &'a HashMap<(u16, u32, u32), BTreeSet<Box<str>>>
    ) -> BTreeSet<&'a str> {
        let mut gene_ids: BTreeSet<&str> = BTreeSet::new();
        for junction in chain.iter() {
            if let Some(hits) = junctions_to_genes.get(&junction.intron_span_key()) {
                gene_ids.extend(hits.iter().map(|gene_id| &**gene_id));
            }
        }
        gene_ids
    }

    // HashMap<junction chain, BTreeSet<gene ID>>
    let chains_to_genes: HashMap<&Vec<SpliceJunction>, BTreeSet<&str>> = chains
        .iter()
        .map(|&chain| (chain, chain_genes(chain, junctions_to_genes)))
        .collect();

    // Step 6. Perform degradation-aware fold.
    // HashMap<junction chain, HashSet<read ID>>
    let mut spliced_clusters: HashMap<Vec<SpliceJunction>, HashSet<usize>> = HashMap::new();
    for &chain in chains.iter() {
        // The single full-length backbone this chain is committed to.
        let backbone: &Vec<SpliceJunction> = if maximal_superchains.contains(chain) {
            // A maximal superchain is its own (only) backbone.
            chain
        } else {
            // Choose the single best maximal superchain containing this truncated chain.
            // Contiguous containment is transitive, so at least one always exists.
            let chain_gene_ids: &BTreeSet<&str> = &chains_to_genes[chain];
            let genes_introduced = |superchain: &Vec<SpliceJunction>| -> usize {
                chains_to_genes[superchain].difference(chain_gene_ids).count()
            };
            chains_to_superchains[chain]
                .iter()
                .copied()
                .filter(|superchain| maximal_superchains.contains(superchain))
                .max_by(|a, b| {
                    genes_introduced(b).cmp(&genes_introduced(a))                                   // fewest genes asserted without evidence
                        .then_with(|| a.len().cmp(&b.len()))                                        // most splice junctions
                        .then_with(|| chains_to_reads[*a].len().cmp(&chains_to_reads[*b].len()))    // most read support
                        .then_with(|| a.cmp(b))                                                              // deterministic tie-break
                })
                .unwrap()
        };

        // Commit this chain's reads to its single backbone.
        let read_ids: &HashSet<usize> = &chains_to_reads[chain];
        spliced_clusters
            .entry(backbone.clone())
            .or_default()
            .extend(read_ids.iter().copied());
    }

    // Step 7. Identify reference transcripts matched to individual reads.
    let mut reference_transcript_ids: HashSet<Box<str>> = HashSet::new();
    for transcript_model in transcript_models.iter() {
        reference_transcript_ids.extend(transcript_model.get_reference_transcript_ids());
    }

    // Step 8. The reference gene loci: the full span of every annotated gene, per chromosome.
    // HashMap<chromosome ID, HashMap<gene ID, (start, end)>>
    let reference_gene_spans: &HashMap<u16, HashMap<Box<str>, (u32, u32)>> = &annotation_index.gene_spans;

    // Step 9. Key every junction-less (unspliced or single-exon) read by the set of annotated gene
    // loci its span overlaps. Reads overlapping no annotated gene are held back for locus
    // clustering instead of collapsing into one bucket.
    let mut unspliced_read_support: HashMap<UnsplicedClusterKey, HashSet<usize>> = HashMap::new();
    // HashMap<chromosome ID, Vec<(read ID, start, end)>>
    let mut unplaced_reads: HashMap<u16, Vec<(usize, u32, u32)>> = HashMap::new();
    for transcript_model in transcript_models.iter() {
        if chained_read_ids.contains(&transcript_model.get_read_id()) {
            continue;   // spliced, or chained on a back-splice alone
        }

        let (chromosome_1, start): (u16, u32) = transcript_model.get_reference_start_position();
        let (chromosome_2, end): (u16, u32) = transcript_model.get_reference_end_position();

        if chromosome_1 != chromosome_2 {
            continue;
        }

        let mut hit_reference_gene_ids: BTreeSet<Box<str>> = BTreeSet::new();
        if let Some(gene_spans) = reference_gene_spans.get(&chromosome_1) {
            for (gene_id, &(gene_start, gene_end)) in gene_spans.iter() {
                if overlaps(gene_start as isize, gene_end as isize, start as isize, end as isize) {
                    hit_reference_gene_ids.insert(gene_id.clone());
                }
            }
        }

        let read_id: usize = transcript_model.get_read_id();

        if hit_reference_gene_ids.is_empty() {
            unplaced_reads
                .entry(chromosome_1)
                .or_default()
                .push((read_id, start, end));
        } else {
            unspliced_read_support
                .entry(UnsplicedClusterKey::Genes(chromosome_1, hit_reference_gene_ids))
                .or_default()
                .insert(read_id);
        }
    }

    // Step 10. Single-linkage locus clustering for reads with no annotated gene.
    for (chromosome_id, mut reads) in unplaced_reads {
        reads.sort_by_key(|&(_, start, end)| (start, end));
        let mut component: usize = 0;
        let mut component_end: u32 = 0;
        for (i, &(read_id, start, end)) in reads.iter().enumerate() {
            if i == 0 {
                component_end = end;
            } else if start > component_end.saturating_add(max_unspliced_locus_gap) {
                component += 1;
                component_end = end;
            } else {
                component_end = component_end.max(end);
            }
            unspliced_read_support
                .entry(UnsplicedClusterKey::Locus(chromosome_id, component))
                .or_default()
                .insert(read_id);
        }
    }

    // Step 11. The map of junctions to (reference gene ID, reference transcript ID) holds
    // every reference transcript. The ones matched to the reads are picked at the lookup.

    // Step 12. Materialize clusters with globally unique IDs.
    let mut splice_junction_clusters: Vec<SpliceJunctionCluster> = Vec::new();
    let mut cluster_id: usize = 0;

    // A cluster below the read floor still materializes when one of its reads
    // carries a DNA-confirmed variant: the floor exists to suppress noise, and
    // DNA evidence is exactly what distinguishes a rare real transcript from it.
    let carries_dna_variant = |read_ids: &HashSet<usize>| -> bool {
        dna_carrier_read_ids.map_or(false, |carriers| !carriers.is_disjoint(read_ids))
    };

    // Spliced RNA clusters.
    // Keys are sorted so cluster IDs are stable across runs, exactly as the
    // unspliced branch below does; `cluster_id` increments inside this loop, so
    // raw HashMap order would renumber every cluster on every run.
    let mut spliced_keys: Vec<&Vec<SpliceJunction>> = spliced_clusters.keys().collect();
    spliced_keys.sort();
    for junctions in spliced_keys {
        let read_ids: &HashSet<usize> = &spliced_clusters[junctions];
        if read_ids.len() < min_reads_per_cluster && !carries_dna_variant(read_ids) {
            continue;
        }

        // Reference gene/transcript junction hits.
        // HashMap<(reference gene ID, reference transcript ID), counter>
        let mut hits: HashMap<(&str, &str), usize> = HashMap::new();
        for junction in junctions.iter() {
            let (chromosome_id, start, end): (u16, u32, u32) = junction.intron_span_key();
            let key: (u16, u32, u32, Strand) = (chromosome_id, start, end, junction.strand_1.clone());
            if let Some(references) = annotation_index.junctions_to_transcripts.get(&key) {
                for &reference in references.iter() {
                    let (reference_gene_id, reference_transcript_id): &(Box<str>, Box<str>) = &annotation_index.transcripts[reference as usize];
                    if reference_transcript_ids.contains(reference_transcript_id) {
                        *hits.entry((reference_gene_id, reference_transcript_id)).or_insert(0) += 1;
                    }
                }
            }
        }

        let mut reference_gene_transcript_ids: Vec<(Box<str>, Box<str>)> =
            match hits.values().copied().max() {
                Some(max_count) => hits
                    .iter()
                    .filter(|&(_, &count)| count == max_count)
                    .map(|(&(gene_id, transcript_id), _)| (gene_id.into(), transcript_id.into()))
                    .collect(),
                None => Vec::new()
            };
        reference_gene_transcript_ids.sort();

        // The chain minus its back-splices: downstream reads this list as introns.
        let splice_junctions: Vec<SpliceJunction> = junctions
            .iter()
            .filter(|junction| !cycle_junctions.contains(*junction))
            .cloned()
            .collect();
        splice_junction_clusters.push(SpliceJunctionCluster::new(
            cluster_id,
            splice_junctions,
            read_ids.clone(),
            reference_gene_transcript_ids
        ));
        cluster_id += 1;
    }

    // Unspliced (junction-less) RNA clusters.
    // Keys are sorted so cluster IDs are stable across runs.
    let mut unspliced_keys: Vec<&UnsplicedClusterKey> = unspliced_read_support.keys().collect();
    unspliced_keys.sort();
    for key in unspliced_keys {
        let read_ids: &HashSet<usize> = &unspliced_read_support[key];
        if read_ids.len() < min_reads_per_cluster && !carries_dna_variant(read_ids) {
            continue;
        }
        let reference_gene_transcript_ids: Vec<(Box<str>, Box<str>)> = match key {
            UnsplicedClusterKey::Genes(_, gene_ids) => gene_ids
                .iter()
                .map(|gene_id| (gene_id.clone(), "".into()))
                .collect(),
            // No annotated gene overlaps this locus, so there is nothing to name it by.
            UnsplicedClusterKey::Locus(_, _) => Vec::new()
        };
        splice_junction_clusters.push(SpliceJunctionCluster::new(
            cluster_id,
            Vec::new(),
            read_ids.clone(),
            reference_gene_transcript_ids
        ));
        cluster_id += 1;
    }

    splice_junction_clusters
}


#[cfg(test)]
#[path = "../tests/clustering/splice_junction_clustering.rs"]
mod tests;