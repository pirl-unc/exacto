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


use edit_distance::edit_distance;
use exacto_core::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::prelude::*;
use crate::calling::dna::variant_record_clustering::calculate_max_dna_clustering_distance;


pub fn cluster_rna_variant_records(
    splice_junctions: &[SpliceJunction],
    transcript_models: &Vec<Arc<TranscriptModel>>,
    max_clustering_distance: u32,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    sequencing_error: f64,
    poa_match_score: i32,
    poa_mismatch_score: i32,
    poa_gap_open_score: i32,
    poa_gap_extend_score: i32
) -> Vec<VariantCall> {
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(splice_junctions);

    // Step 1. Route variant records by type.
    let mut breakpoint_variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    let mut nonbreakpoint_variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    for variant_record in transcript_models.iter().flat_map(|transcript_model| transcript_model.get_variant_records().iter()) {
        match variant_record.get_variant_type() {
            VariantType::Breakpoint
            | VariantType::CircularRNA
            | VariantType::FusionGene
            | VariantType::Translocation => {
                breakpoint_variant_records.push(Arc::new(variant_record.clone()))
            },
            VariantType::SingleNucleotideVariant
            | VariantType::MultiNucleotideVariant
            | VariantType::Insertion
            | VariantType::Deletion=> {
                nonbreakpoint_variant_records.push(Arc::new(variant_record.clone()))
            },
            VariantType::CrypticExon
            | VariantType::ExonTruncation
            | VariantType::IntronRetention
            | VariantType::NonCanonicalSplicing
            | VariantType::UTRExtension => {
                // Do not cluster
            }
        }
    }

    // Step 2. Cluster the breakpoint variant records.
    let mut variant_record_clusters: Vec<VariantRecordCluster> = cluster_breakpoint_rna_variant_records(
        breakpoint_variant_records,
        &ladders,
        max_clustering_distance
    );

    // Step 3. Cluster the non-breakpoint variant records.
    let mut variant_records_by_chromosome: HashMap<(ReferenceChromosomeID, ReferenceChromosomeID), Vec<Arc<VariantRecord>>> = HashMap::new();
    for record in nonbreakpoint_variant_records {
        variant_records_by_chromosome
            .entry((record.get_chromosome_1(), record.get_chromosome_2()))
            .or_default()
            .push(record);
    }
    for records in variant_records_by_chromosome.values() {
        variant_record_clusters.extend(cluster_non_breakpoint_rna_variant_records(
            records,
            &ladders,
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_clustering_distance,
            sequencing_error
        ));
    }

    // Step 4. Create variant calls.
    let mut variant_calls: Vec<VariantCall> = Vec::new();
    for variant_record_cluster in variant_record_clusters.iter() {
        // Get the variant records from the cluster.
        let variant_records: HashSet<VariantRecord> = variant_record_cluster
            .get_variant_records()
            .iter()
            .map(Arc::as_ref)
            .cloned()
            .collect();

        // Create a variant call.
        let variant_call: VariantCall = VariantCall::from_variant_records(
            0,
            variant_records,
            poa_match_score,
            poa_mismatch_score,
            poa_gap_open_score,
            poa_gap_extend_score
        );

        variant_calls.push(variant_call);
    }

    variant_calls
}


/// Which breakpoint of a two-sided record a pool is built on.
#[derive(Clone, Copy)]
enum Side {
    One,
    Two
}


pub(crate) struct ExonLadder {
    intron_starts: Vec<ReferencePosition>,
    intron_ends: Vec<ReferencePosition>
}

impl ExonLadder {
    pub fn from_splice_junctions(splice_junctions: &[SpliceJunction], chromosome: ReferenceChromosomeID) -> Self {
        let mut introns: Vec<(ReferencePosition, ReferencePosition)> = splice_junctions
            .iter()
            .filter(|junction| junction.chromosome_1 == chromosome && junction.chromosome_2 == chromosome)
            .map(SpliceJunction::intron_span)
            .collect();
        introns.sort_unstable();
        introns.dedup();
        Self {
            intron_starts: introns.iter().map(|i| i.0).collect(),
            intron_ends: introns.iter().map(|i| i.1).collect()
        }
    }

    pub fn exon_index(&self, position: ReferencePosition) -> usize {
        self.intron_starts.partition_point(|&start| start < position)
    }
    
    pub fn covers(&self, position: ReferencePosition) -> bool {
        match (self.intron_starts.first(), self.intron_ends.last()) {
            (Some(&first), Some(&last)) => position >= first && position <= last,
            _ => false
        }
    }

    pub fn exon_distance(&self, a: ReferencePosition, b: ReferencePosition) -> usize {
        self.exon_index(a).abs_diff(self.exon_index(b))
    }

    pub fn intron_end(&self, position: ReferencePosition) -> Option<ReferencePosition> {
        let k: usize = self.intron_ends.partition_point(|&end| end < position);
        (k < self.intron_starts.len() && self.intron_starts[k] <= position).then(|| self.intron_ends[k])
    }
    
    pub fn spliced_distance(&self, a: ReferencePosition, b: ReferencePosition) -> u32 {
        let (low, high): (ReferencePosition, ReferencePosition) = (a.min(b), a.max(b));
        if low == high {
            return 0;
        }
        
        let first: usize = self.intron_starts.partition_point(|&start| start <= low);
        let last: usize = self.intron_ends.partition_point(|&end| end < high);
        let mut intronic: u32 = 0;
        for k in first..last {
            intronic += self.intron_ends[k] - self.intron_starts[k] + 1;
        }
        (high - low) - intronic
    }
}


pub(crate) struct ExonLadders {
    by_chromosome: HashMap<ReferenceChromosomeID, ExonLadder>
}

impl ExonLadders {
    pub(crate) fn from_splice_junctions(splice_junctions: &[SpliceJunction]) -> Self {
        let chromosomes: HashSet<ReferenceChromosomeID> = splice_junctions
            .iter()
            .flat_map(|junction| [junction.chromosome_1, junction.chromosome_2])
            .collect();
        Self {
            by_chromosome: chromosomes
                .into_iter()
                .map(|chromosome| (chromosome, ExonLadder::from_splice_junctions(splice_junctions, chromosome)))
                .collect()
        }
    }
    
    pub fn exon_distance(&self, chromosome: ReferenceChromosomeID, a: ReferencePosition, b: ReferencePosition) -> Option<usize> {
        let ladder: &ExonLadder = self.by_chromosome.get(&chromosome)?;
        (ladder.covers(a) && ladder.covers(b)).then(|| ladder.exon_distance(a, b))
    }
    
    pub fn spliced_distance(&self, chromosome: ReferenceChromosomeID, a: ReferencePosition, b: ReferencePosition) -> u32 {
        match self.by_chromosome.get(&chromosome) {
            Some(ladder) => ladder.spliced_distance(a, b),
            None => a.abs_diff(b)
        }
    }
    
    pub fn next_sweep_index(&self, chromosome: ReferenceChromosomeID, a: ReferencePosition, positions: &[ReferencePosition], index: usize) -> Option<usize> {
        let ladder: &ExonLadder = self.by_chromosome.get(&chromosome)?;
        let intron_end: ReferencePosition = ladder.intron_end(positions[index])?;
        if ladder.intron_end(a) == Some(intron_end) {
            return None;
        }
        Some(index + positions[index..].partition_point(|&position| position <= intron_end))
    }
}


/// (chromosome, position, operation type) of one side of a record.
fn breakend(record: &VariantRecord, side: Side) -> (ReferenceChromosomeID, ReferencePosition, GraphOperationType) {
    match side {
        Side::One => (record.get_chromosome_1(), record.get_position_1(), record.get_operation_1().clone()),
        Side::Two => (record.get_chromosome_2(), record.get_position_2(), record.get_operation_2().clone())
    }
}


pub fn breakpoint_label_rank(variant_type: &VariantType) -> u8 {
    match variant_type {
        VariantType::FusionGene => 3,
        VariantType::CircularRNA => 2,
        VariantType::Translocation => 1,
        VariantType::Breakpoint => 0,
        other => panic!("{other:?} is not a breakpoint-family variant type.")
    }
}


fn breakend_position(record: &VariantRecord, side: Side) -> ReferencePosition {
    match side {
        Side::One => record.get_position_1(),
        Side::Two => record.get_position_2()
    }
}


fn cluster_breakpoint_rna_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    ladders: &ExonLadders,
    max_breakend_distance: u32
) -> Vec<VariantRecordCluster> {
    assert!(
        variant_records.iter().all(|record| matches!(
            record.get_variant_type(),
            VariantType::Breakpoint
            | VariantType::CircularRNA
            | VariantType::FusionGene
            | VariantType::Translocation
        )),
        "Only breakpoint variant types (breakpoint, circular RNA, fusion gene, and translocation) are expected here."
    );

    assert!(
        variant_records.iter().all(|record| *record.get_operation_2() != GraphOperationType::Noop),
        "RNA breakpoint records are always resolved."
    );

    if variant_records.is_empty() {
        return Vec::new();
    }

    // Step 1. Deterministic input order; everything downstream keys on indices.
    let mut variant_records: Vec<Arc<VariantRecord>> = variant_records;
    variant_records.sort_by_cached_key(|record| (
        record.get_read_id(),
        record.get_read_position_1(),
        record.get_read_position_2(),
        record.get_graph_operation().as_boxed_str()
    ));

    // Step 2. Pool each side independently, in spliced coordinates.
    let pool_1: Vec<usize> = pool_breakpoints(
        &variant_records,
        Side::One,
        ladders,
        max_breakend_distance
    );
    let pool_2: Vec<usize> = pool_breakpoints(
        &variant_records,
        Side::Two,
        ladders,
        max_breakend_distance
    );

    // Step 3. A junction is a pair of pools. BTreeMap keeps the output order deterministic.
    let mut junctions: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for i in 0..variant_records.len() {
        junctions.entry((pool_1[i], pool_2[i])).or_default().push(i);
    }

    // Step 4. One cluster per junction, members in input order, all under the junction's
    // most specific label.
    junctions
        .into_values()
        .map(|indices| {
            let variant_type: VariantType = indices
                .iter()
                .map(|&i| variant_records[i].get_variant_type())
                .max_by_key(|variant_type| breakpoint_label_rank(variant_type))
                .unwrap()
                .clone();
            let members: Vec<Arc<VariantRecord>> = indices
                .iter()
                .map(|&i| {
                    let record: &Arc<VariantRecord> = &variant_records[i];
                    if *record.get_variant_type() == variant_type {
                        Arc::clone(record)
                    } else {
                        Arc::new(VariantRecord::new(
                            record.get_read_id(),
                            record.get_read_position_1(),
                            record.get_read_position_2(),
                            record.get_graph_operation().with_variant_type(variant_type.clone())
                        ))
                    }
                })
                .collect();
            VariantRecordCluster::from_variant_records(&members)
        })
        .collect()
}


fn cluster_non_breakpoint_rna_variant_records(
    variant_records: &Vec<Arc<VariantRecord>>,
    ladders: &ExonLadders,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64
) -> Vec<VariantRecordCluster> {
    assert!(
        variant_records.iter().all(|record| !matches!(
            record.get_variant_type(),
            VariantType::Breakpoint | VariantType::Translocation
        ))
    );

    if variant_records.is_empty() {
        return Vec::new();
    }

    // Step 1. One chromosome pair per call; the caller splits on it.
    let chromosome_1: ReferenceChromosomeID = variant_records[0].get_chromosome_1();
    let chromosome_2: ReferenceChromosomeID = variant_records[0].get_chromosome_2();
    assert!(
        variant_records.
            iter()
            .all(|r| r.get_chromosome_1() == chromosome_1 && r.get_chromosome_2() == chromosome_2)
    );

    // Step 2. Genomic order on position 1; spliced distance is monotone in it along the exons.
    let mut variant_records: Vec<Arc<VariantRecord>> = variant_records.clone();
    variant_records.sort_by_key(|r| (r.get_position_1(), r.get_position_2(), r.get_read_id()));

    // Step 3. Group variant records by Grammar and Spelling. 
    type Grammar = (VariantType, GraphOperationType, GraphOperationType, ReferenceChromosomeID, ReferenceChromosomeID, Strand, Strand);
    // (reference position, reference position, sequence)
    type Spelling = (ReferencePosition, ReferencePosition, Box<str>);
    // HashMap<Grammar, HashMap<Spelling, Vec<variant record index>>>
    let mut groups: HashMap<Grammar, HashMap<Spelling, Vec<usize>>> = HashMap::new();
    for (i, record) in variant_records.iter().enumerate() {
        let (strand_1, strand_2): (Strand, Strand) = grammar_strands(record);
        let grammar: Grammar = (
            record.get_variant_type().clone(),
            record.get_operation_1().clone(),
            record.get_operation_2().clone(),
            record.get_chromosome_1(),
            record.get_chromosome_2(),
            strand_1,
            strand_2
        );
        let spelling: Spelling = (
            record.get_position_1(),
            record.get_position_2(),
            record.get_standardized_sequence().into()
        );
        groups
            .entry(grammar)
            .or_insert_with(HashMap::new)
            .entry(spelling)
            .or_insert_with(Vec::new)
            .push(i);
    }

    // Step 4. Representatives.
    let mut representatives: Vec<Arc<VariantRecord>> = Vec::new();
    let mut members_of: Vec<Vec<usize>> = Vec::new();
    let mut partitions: Vec<Vec<usize>> = Vec::new();
    let mut synthetic_read_id: usize = usize::MAX;
    for spellings in groups.into_values() {
        let mut partition: Vec<usize> = Vec::new();
        for members in spellings.into_values() {
            let first: &Arc<VariantRecord> = &variant_records[members[0]];
            let one_read: bool = members
                .iter()
                .all(|&i| variant_records[i].get_read_id() == first.get_read_id());
            if one_read {
                for &i in members.iter() {
                    partition.push(representatives.len());
                    representatives.push(Arc::clone(&variant_records[i]));
                    members_of.push(vec![i]);
                }
            } else {
                partition.push(representatives.len());
                representatives.push(Arc::new(VariantRecord::new(
                    synthetic_read_id,
                    first.get_read_position_1(),
                    first.get_read_position_2(),
                    first.get_graph_operation().clone()
                )));
                synthetic_read_id -= 1;
                members_of.push(members);
            }
        }
        partitions.push(partition);
    }

    // Step 5. The sweep, per grammar, over representatives: same sort key, same criteria,
    // break at the grammar's own tolerance.
    let mut uf: UnionFind = UnionFind::new();
    for partition in partitions.iter_mut() {
        partition.sort_by_key(|&r| {
            let record: &VariantRecord = representatives[r].as_ref();
            (record.get_position_1(), record.get_position_2(), record.get_read_id())
        });
        let sizes: Vec<isize> = partition
            .iter()
            .map(|&r| representatives[r].get_variant_size())
            .collect();
        let window: u32 = if sizes.iter().all(|&size| size >= 0) {
            match *sizes.iter().max().unwrap() as u32 {
                1 => 0,
                size => calculate_max_dna_clustering_distance(
                    size,
                    sequencing_error,
                    max_clustering_distance
                )
            }
        } else {
            max_clustering_distance
        };
        let positions: Vec<ReferencePosition> = partition
            .iter()
            .map(|&r| representatives[r].get_position_1())
            .collect();
        for (k, &ri) in partition.iter().enumerate() {
            uf.union(ri as u32, ri as u32);
            let a: &VariantRecord = representatives[ri].as_ref();
            let mut m: usize = k + 1;
            while m < partition.len() {
                let rj: usize = partition[m];
                let b: &VariantRecord = representatives[rj].as_ref();
                if ladders.spliced_distance(chromosome_1, a.get_position_1(), b.get_position_1()) > window {
                    match ladders.next_sweep_index(chromosome_1, a.get_position_1(), &positions, m) {
                        Some(next) => {
                            m = next;
                            continue;
                        }
                        None => break
                    }
                }
                if meet_rna_clustering_criteria(
                    a,
                    b,
                    ladders,
                    min_size_proportion,
                    max_ins_norm_edit_distance,
                    max_clustering_distance,
                    sequencing_error
                ) {
                    uf.union(ri as u32, rj as u32);
                }
                m += 1;
            }
        }
    }

    // Step 6. Expand representatives to their records; clusters ordered by their first
    // member, as before.
    let mut groups: Vec<Vec<usize>> = uf
        .get_clusters()
        .into_iter()
        .map(|g| {
            let mut g: Vec<usize> = g
                .into_iter()
                .flat_map(|r| members_of[r as usize].iter().copied())
                .collect();
            g.sort_unstable();
            g
        })
        .collect();
    groups.sort_unstable_by_key(|g| g[0]);

    groups
        .into_iter()
        .map(|indices| {
            let members: Vec<Arc<VariantRecord>> = indices
                .iter()
                .map(|&i| Arc::clone(&variant_records[i]))
                .collect();
            VariantRecordCluster::from_variant_records(&members)
        })
        .collect()
}


fn grammar_strands(record: &VariantRecord) -> (Strand, Strand) {
    match record.get_variant_type() {
        VariantType::SingleNucleotideVariant
        | VariantType::MultiNucleotideVariant
        | VariantType::Insertion
        | VariantType::Deletion => (Strand::Unknown, Strand::Unknown),
        _ => (record.get_strand_1().clone(), record.get_strand_2().clone())
    }
}


fn meet_rna_clustering_criteria(
    a: &VariantRecord,
    b: &VariantRecord,
    ladders: &ExonLadders,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64
) -> bool {
    // One read never supports a site twice.
    if a.get_read_id() == b.get_read_id() {
        return false;
    }
    // Same grammar on both sides.
    if a.get_variant_type() != b.get_variant_type()
        || a.get_operation_1() != b.get_operation_1()
        || a.get_operation_2() != b.get_operation_2()
        || a.get_chromosome_1() != b.get_chromosome_1()
        || a.get_chromosome_2() != b.get_chromosome_2()
        || grammar_strands(a) != grammar_strands(b) {
        return false;
    }

    // Size proportion, when both sizes are meaningful.
    let a_size: isize = a.get_variant_size();
    let b_size: isize = b.get_variant_size();
    let has_valid_sizes: bool = a_size >= 0 && b_size >= 0;
    if has_valid_sizes {
        let max_size: f64 = a_size.max(b_size) as f64;
        if max_size > 0.0 && (a_size.min(b_size) as f64) / max_size < min_size_proportion {
            return false;
        }
    }

    // Proximity in spliced coordinates, tolerance scaled by size as for DNA.
    let max_distance: u32 = if has_valid_sizes {
        match a_size.max(b_size) as u32 {
            1 => 0,
            size => calculate_max_dna_clustering_distance(size, sequencing_error, max_clustering_distance)
        }
    } else {
        max_clustering_distance
    };
    let distance_1: u32 = ladders.spliced_distance(a.get_chromosome_1(), a.get_position_1(), b.get_position_1());
    let distance_2: u32 = ladders.spliced_distance(a.get_chromosome_2(), a.get_position_2(), b.get_position_2());
    if distance_1 > max_distance || distance_2 > max_distance {
        return false;
    }

    // Insertions must also agree on sequence.
    if *a.get_variant_type() == VariantType::Insertion {
        let a_sequence: String = a.get_standardized_sequence();
        let b_sequence: String = b.get_standardized_sequence();
        let max_length: f64 = a_sequence.len().max(b_sequence.len()) as f64;
        let normalized: f64 = edit_distance(&a_sequence, &b_sequence) as f64 / max_length;
        if normalized > max_ins_norm_edit_distance {
            return false;
        }
    }

    // Substitutions must spell the same bases: another base at one position is another allele.
    if matches!(a.get_variant_type(), VariantType::SingleNucleotideVariant | VariantType::MultiNucleotideVariant)
        && a.get_standardized_sequence() != b.get_standardized_sequence() {
        return false;
    }
    true
}


fn pool_breakpoints(
    variant_records: &[Arc<VariantRecord>],
    side: Side,
    ladders: &ExonLadders,
    max_breakend_distance: u32
) -> Vec<usize> {
    // Step 1. Group record indices by (chromosome, operation type).
    let mut groups: HashMap<(ReferenceChromosomeID, GraphOperationType), Vec<usize>> = HashMap::new();
    for (i, record) in variant_records.iter().enumerate() {
        let (chromosome, _, operation) = breakend(record, side);
        groups.entry((chromosome, operation)).or_default().push(i);
    }

    // Step 2. Sweep each group by position.
    let mut pools: Vec<Vec<usize>> = Vec::new();
    for ((chromosome, _), indices) in groups.iter_mut() {
        indices.sort_unstable_by_key(|&i| breakend_position(&variant_records[i], side));
        let positions: Vec<ReferencePosition> = indices
            .iter()
            .map(|&i| breakend_position(&variant_records[i], side))
            .collect();
        let mut uf: UnionFind = UnionFind::new();
        for (k, &i) in indices.iter().enumerate() {
            uf.union(i as u32, i as u32);
            let a: ReferencePosition = positions[k];
            let mut m: usize = k + 1;
            while m < indices.len() {
                if ladders.spliced_distance(*chromosome, a, positions[m]) > max_breakend_distance {
                    match ladders.next_sweep_index(*chromosome, a, &positions, m) {
                        Some(next) => {
                            m = next;
                            continue;
                        }
                        None => break
                    }
                }
                uf.union(i as u32, indices[m] as u32);
                m += 1;
            }
        }
        pools.extend(
            uf.get_clusters()
                .into_iter()
                .map(|pool| pool.into_iter().map(|i| i as usize).collect::<Vec<usize>>())
        );
    }

    // Step 3. Deterministic pool ids: order pools by their smallest member.
    pools.sort_unstable_by_key(|pool| *pool.iter().min().unwrap());
    let mut pool_of: Vec<usize> = vec![0; variant_records.len()];
    for (pool_id, pool) in pools.iter().enumerate() {
        for &i in pool.iter() {
            pool_of[i] = pool_id;
        }
    }
    pool_of
}


#[cfg(test)]
#[path = "../../tests/calling/rna/variant_record_clustering.rs"]
mod tests;