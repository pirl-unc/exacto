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
use exacto_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::cmp::Reverse;

use crate::prelude::*;


pub fn identify_reference_transcript_matches(
    query_exons: &Vec<TranscriptModelExon>,
    query_splice_junctions: &Vec<TranscriptModelSpliceJunction>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>
) -> Vec<ReferenceTranscriptMatch> {
    // Step 1. Identifying overlapping reference transcript IDs.
    let reference_transcript_ids: HashSet<ReferenceTranscriptID> = identify_overlapping_reference_transcript_ids(
        query_exons,
        gene_annotator,
        chromosome_names_map
    );

    // Step 2. Identifying overlapping reference gene IDs.
    let mut reference_gene_ids: HashSet<ReferenceGeneID> = HashSet::new();
    for reference_transcript_id in reference_transcript_ids.iter() {
        let transcript: &Transcript = gene_annotator.get_transcript(reference_transcript_id).unwrap();
        reference_gene_ids.insert(transcript.gene_id.clone());
    }

    // Step 3. Score each reference transcript in each reference gene.
    let mut reference_transcript_scores: HashMap<ReferenceGeneID, Vec<ReferenceTranscriptMatch>> = HashMap::new();
    for reference_gene_id in reference_gene_ids.iter() {
        let reference_gene: &Gene = gene_annotator.get_gene(reference_gene_id).unwrap();
        let mut reference_transcripts: Vec<&Transcript> = Vec::new();
        for reference_transcript in reference_gene.transcripts.values() {
            if reference_transcript_ids.contains(&reference_transcript.transcript_id) {
                reference_transcripts.push(reference_transcript);
            }
        }
        let reference_transcript_matches_: Vec<ReferenceTranscriptMatch> = score_reference_transcripts(
            query_exons,
            query_splice_junctions,
            reference_transcripts,
            reference_gene,
            chromosome_names_map
        );
        reference_transcript_scores.insert(reference_gene.gene_id.clone(), reference_transcript_matches_);
    }

    // Step 4. Select the best reference transcript match for each gene.
    let mut reference_transcript_matches: Vec<ReferenceTranscriptMatch> = Vec::new();
    for reference_gene_id in reference_gene_ids.iter() {
        if let Some(matches) = reference_transcript_scores.get(reference_gene_id) {
            // Build reference transcript canonical ordering.
            let ranked: Vec<&Transcript> = gene_annotator.rank_transcripts(
                matches.iter()
                    .map(|m| gene_annotator.get_transcript(&*m.get_reference_transcript_id()).unwrap())
                    .collect()
            );
            let canonical_order: HashMap<&str, usize> = ranked.iter()
                .enumerate()
                .map(|(i, t)| (t.transcript_id.as_ref(), i))
                .collect();

            // Select the best reference transcript match.
            let best_match: Option<&ReferenceTranscriptMatch> = matches
                .iter()
                .min_by_key(|m| rank_reference_transcript_match(m, gene_annotator, &canonical_order));

            assert!(best_match.is_some());

            reference_transcript_matches.push(best_match.unwrap().clone());
        }
    }

    // Step 5. Order the matches best first.
    let ranked: Vec<&Transcript> = gene_annotator.rank_transcripts(
        reference_transcript_matches
            .iter()
            .map(|m| gene_annotator.get_transcript(m.get_reference_transcript_id()).unwrap())
            .collect()
    );
    let canonical_order: HashMap<&str, usize> = ranked
        .iter()
        .enumerate()
        .map(|(i, t)| (t.transcript_id.as_ref(), i))
        .collect();
    reference_transcript_matches.sort_by(|a, b| {
        rank_reference_transcript_match(a, gene_annotator, &canonical_order)
            .cmp(&rank_reference_transcript_match(b, gene_annotator, &canonical_order))
    });

    reference_transcript_matches
}


pub fn score_reference_transcript(
    query_exons: &Vec<TranscriptModelExon>,
    query_splice_junctions: &Vec<TranscriptModelSpliceJunction>,
    reference_transcript: &Transcript,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>
) -> ReferenceTranscriptMatch {
    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    struct SpliceSite {
        chromosome: ReferenceChromosomeName,
        position_1: ReferencePosition,
        position_2: ReferencePosition,
        strand: Strand,
    }
    impl SpliceSite {
        fn new(chromosome: &str, position_1: ReferencePosition, position_2: ReferencePosition, strand: &Strand) -> Self {
            Self { chromosome: chromosome.into(), position_1, position_2, strand: strand.clone() }
        }
    }

    // Step 1. Get the assembled transcript's splice junctions.
    let mut query_splice_sites: HashSet<SpliceSite> = HashSet::new();
    for query_splice_junction in query_splice_junctions.iter() {
        if query_splice_junction.reference_chromosome_id_1 == query_splice_junction.reference_chromosome_id_2 
            && query_splice_junction.reference_strand_1 == query_splice_junction.reference_strand_2 {
            let chromosome: ReferenceChromosomeName = chromosome_names_map
                .get_by_right(&query_splice_junction.reference_chromosome_id_1)
                .unwrap()
                .clone();
            let (position_1, position_2) = if query_splice_junction.reference_position_1 < query_splice_junction.reference_position_2 {
                (query_splice_junction.reference_position_1, query_splice_junction.reference_position_2)
            } else {
                (query_splice_junction.reference_position_2, query_splice_junction.reference_position_1)
            };
            let sj: SpliceSite = SpliceSite::new(
                &*chromosome,
                position_1,
                position_2,
                &query_splice_junction.reference_strand_1
            );
            query_splice_sites.insert(sj);
        }
    }


    // Step 2. Get the reference transcript's splice junctions (i.e. introns).
    let reference_introns: Vec<Intron> = reference_transcript.get_introns();
    let mut reference_splice_sites: HashSet<SpliceSite> = HashSet::new();
    for reference_intron in reference_introns.iter() {
        let sj: SpliceSite = SpliceSite::new(
            &*reference_intron.chromosome,
            reference_intron.start,
            reference_intron.end,
            &reference_intron.strand
        );
        reference_splice_sites.insert(sj);
    }


    // Step 3. Identify the matched splice junctions.
    let splice_junction_matches: Vec<ReferenceTranscriptSpliceJunction> = reference_splice_sites
        .intersection(&query_splice_sites)
        .map(|s| {
            ReferenceTranscriptSpliceJunction::new(
                &*s.chromosome,
                s.position_1,
                s.position_2,
                s.strand.clone()
            )
        })
        .collect();

    // Step 4. Get the transcript model's exonic regions.
    let model_exon_regions: Vec<(ReferenceChromosomeName, ReferencePosition, ReferencePosition)> = query_exons
        .iter()
        .map(|exon| {
            let chr = chromosome_names_map
                .get_by_right(&exon.reference_chromosome_id)
                .unwrap()
                .to_string()
                .into_boxed_str();
            (chr, exon.reference_start, exon.reference_end)
        })
        .collect();


    // Step 5. Get the reference transcript's exonic regions.
    let reference_exons: Vec<(ReferenceChromosomeName, ReferencePosition, ReferencePosition)> = reference_transcript
        .get_sorted_exons()
        .iter()
        .map(|exon| (exon.chromosome.clone(), exon.start, exon.end))
        .collect();


    // Step 6. Count the shared and private bases.
    let num_overlapping_bases: u32 = count_common_bases(&model_exon_regions, &reference_exons);
    let (num_query_only_bases, num_reference_only_bases) = count_non_overlapping_bases(
        &model_exon_regions,
        &reference_exons
    );


    // Step 7. Split the query-only bases at the reference transcript's span.
    let reference_span: Vec<(ReferenceChromosomeName, ReferencePosition, ReferencePosition)> = match reference_exons.first() {
        Some((chromosome, _, _)) => vec![(
            chromosome.clone(),
            reference_exons.iter().map(|(_, start, _)| *start).min().unwrap(),
            reference_exons.iter().map(|(_, _, end)| *end).max().unwrap()
        )],
        None => Vec::new()
    };
    let (num_terminal_query_only_bases, _) = count_non_overlapping_bases(
        &model_exon_regions,
        &reference_span
    );


    let reference_transcript_match: ReferenceTranscriptMatch = ReferenceTranscriptMatch::new(
        &*reference_gene.gene_id,
        &*reference_gene.gene_name,
        &*reference_transcript.transcript_id,
        &splice_junction_matches,
        num_overlapping_bases,
        num_reference_only_bases,
        num_query_only_bases,
        num_terminal_query_only_bases
    );

    reference_transcript_match
}


/// Rank a reference transcript match. The smaller rank is the better match:
///  1) Highest splice junction match count, then
///  2) Smallest number of TERMINAL query-only bases, then
///  3) Smallest number of INTERNAL query-only bases, then
///  4) protein_coding before every other biotype, then
///  5) Smallest number of reference-only bases, then
///  6) Canonical ordering (lower rank wins), then
///  7) Reference transcript ID, so the order is total.
fn rank_reference_transcript_match<'a>(
    reference_transcript_match: &'a ReferenceTranscriptMatch,
    gene_annotator: &impl GeneAnnotator,
    canonical_order: &HashMap<&str, usize>
) -> (Reverse<usize>, u32, u32, bool, u32, usize, &'a str) {
    let is_protein_coding: bool = gene_annotator
        .get_transcript(reference_transcript_match.get_reference_transcript_id())
        .map(|t| &*t.transcript_type == "protein_coding")
        .unwrap_or(false);
    (
        Reverse(reference_transcript_match.num_splice_junction_matches()),                                      // higher count wins
        reference_transcript_match.num_terminal_query_only_bases(),                                             // fewer past the annotation's ends
        reference_transcript_match.num_internal_query_only_bases(),                                             // fewer retained intron bases
        !is_protein_coding,                                                                                     // protein_coding first
        reference_transcript_match.num_reference_only_bases(),                                                  // fewer reference-only
        *canonical_order.get(reference_transcript_match.get_reference_transcript_id()).unwrap_or(&usize::MAX),  // canonical rank
        reference_transcript_match.get_reference_transcript_id()                                                // total order: no tie left to iteration order
    )
}


fn identify_overlapping_reference_transcript_ids(
    model_exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>
) -> HashSet<ReferenceTranscriptID> {
    let mut reference_transcript_ids: HashSet<ReferenceTranscriptID> = HashSet::new();
    for model_exon in model_exons.iter() {
        let chromosome_name: ReferenceChromosomeName = chromosome_names_map
            .get_by_right(&model_exon.reference_chromosome_id)
            .unwrap()
            .to_string()
            .into_boxed_str();
        let overlapping_transcript_ids: Vec<ReferenceTranscriptID> = gene_annotator.get_transcript_ids_overlapping_region(
            &*chromosome_name,
            model_exon.reference_start,
            model_exon.reference_end
        );
        for transcript_id in overlapping_transcript_ids.iter() {
            if reference_transcript_ids.contains(transcript_id) {
                continue;
            }
            let transcript: &Transcript = gene_annotator.get_transcript(transcript_id).unwrap();
            if transcript.strand == model_exon.reference_strand {
                for reference_exon in transcript.get_sorted_exons() {
                    if overlaps(
                        model_exon.reference_start as isize,
                        model_exon.reference_end as isize,
                        reference_exon.start as isize,
                        reference_exon.end as isize
                    ) {
                        reference_transcript_ids.insert(transcript_id.clone());
                    }
                }
            }
        }
    }

    reference_transcript_ids
}


fn score_reference_transcripts(
    query_exons: &Vec<TranscriptModelExon>,
    query_splice_junctions: &Vec<TranscriptModelSpliceJunction>,
    reference_transcripts: Vec<&Transcript>,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>
) -> Vec<ReferenceTranscriptMatch> {
    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = reference_transcripts
        .iter()
        .map(|reference_transcript| {
            let reference_transcript_match: ReferenceTranscriptMatch = score_reference_transcript(
                query_exons,
                query_splice_junctions,
                reference_transcript,
                reference_gene,
                chromosome_names_map
            );
            reference_transcript_match
        })
        .collect::<Vec<_>>();
    reference_transcript_matches
}


#[cfg(test)]
#[path = "../tests/reference/reference_transcript_matching.rs"]
mod tests;