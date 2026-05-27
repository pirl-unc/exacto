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

use crate::prelude::*;


pub fn identify_reference_transcript_matches(
    exons: &Vec<AssembledTranscriptExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    scoring_method: ReferenceTranscriptScoringMethod,
    selection_strategy: ReferenceTranscriptSelectionStrategy,
    top_k: usize,
    threshold: f32
) -> Vec<ReferenceTranscriptMatch> {
    // Step 1. Identifying overlapping transcript IDa
    let reference_transcript_ids: HashSet<Box<str>> = identify_overlapping_transcript_ids(
        exons,
        gene_annotator,
        chromosome_names_map
    );

    // Step 2. Identifying overlapping gene IDs
    let mut reference_gene_ids: HashSet<Box<str>> = HashSet::new();
    for reference_transcript_id in reference_transcript_ids.iter() {
        let transcript: &Transcript = gene_annotator.get_transcript(reference_transcript_id).unwrap();
        reference_gene_ids.insert(transcript.gene_id.clone());
    }

    // Step 3. Score each overlapping transcript in each gene
    let mut reference_transcript_scores: HashMap<Box<str>, Vec<ReferenceTranscriptMatch>> = HashMap::new();
    for reference_gene_id in reference_gene_ids.iter() {
        let reference_gene: &Gene = gene_annotator.get_gene(reference_gene_id).unwrap();
        let mut reference_transcripts: Vec<&Transcript> = Vec::new();
        for reference_transcript in reference_gene.transcripts.values() {
            if reference_transcript_ids.contains(&reference_transcript.transcript_id) {
                reference_transcripts.push(reference_transcript);
            }
        }
        let reference_transcript_matches_: Vec<ReferenceTranscriptMatch> = score_reference_transcripts(
            exons,
            reference_transcripts,
            reference_gene,
            chromosome_names_map,
            scoring_method.clone()
        );
        reference_transcript_scores.insert(reference_gene.gene_id.clone(), reference_transcript_matches_);
    }

    // Step 4. Select reference transcript matches
    let mut reference_transcript_matches: Vec<ReferenceTranscriptMatch> = Vec::new();
    for reference_gene_id in reference_gene_ids.iter() {
        if let Some(matches) = reference_transcript_scores.get(reference_gene_id) {
            let mut sorted_matches: Vec<ReferenceTranscriptMatch> = matches.clone();

            // Build canonical ordering via rank_transcripts
            let ranked: Vec<&Transcript> = gene_annotator.rank_transcripts(
                sorted_matches.iter()
                    .map(|m| gene_annotator.get_transcript(&*m.reference_transcript_id).unwrap())
                    .collect()
            );
            let canonical_order: HashMap<&str, usize> = ranked.iter()
                .enumerate()
                .map(|(i, t)| (t.transcript_id.as_ref(), i))
                .collect();

            // Sort by score descending, then canonical rank ascending as tiebreaker
            sorted_matches.sort_by(|a, b| {
                b.score.partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        let rank_a = canonical_order.get(a.reference_transcript_id.as_ref()).unwrap_or(&usize::MAX);
                        let rank_b = canonical_order.get(b.reference_transcript_id.as_ref()).unwrap_or(&usize::MAX);
                        rank_a.cmp(rank_b)
                    })
            });

            if selection_strategy == ReferenceTranscriptSelectionStrategy::TopK {
                sorted_matches.truncate(top_k);
                reference_transcript_matches.extend(sorted_matches);
            } else if selection_strategy == ReferenceTranscriptSelectionStrategy::Threshold {
                reference_transcript_matches.extend(
                    sorted_matches.into_iter().filter(|m| m.score >= threshold)
                );
            } else {
                panic!("Unsupported selection strategy: {}", selection_strategy.as_str());
            }
        }
    }

    // Step 5. Sort by score
    reference_transcript_matches.sort_by(|a, b|  {
        b.score
            .partial_cmp(&a.score)
            .unwrap()
            .then_with(|| a.reference_transcript_id.cmp(&b.reference_transcript_id))
    });

    reference_transcript_matches
}

fn identify_overlapping_transcript_ids(
    model_exons: &Vec<AssembledTranscriptExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>, u16>
) -> HashSet<Box<str>> {
    let mut reference_transcript_ids: HashSet<Box<str>> = HashSet::new();
    for model_exon in model_exons.iter() {
        let chromosome_name: Box<str> = chromosome_names_map.get_by_right(&model_exon.reference_chromosome_id).unwrap().to_string().into_boxed_str();
        let overlapping_transcript_ids: Vec<Box<str>> = gene_annotator.get_transcript_ids_overlapping_region(
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

fn score_reference_transcript(
    exons: &Vec<AssembledTranscriptExon>,
    reference_transcript: &Transcript,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    scoring_method: ReferenceTranscriptScoringMethod
) -> ReferenceTranscriptMatch {
    // Step 1. Get the transcript model's exonic regions
    let model_exon_regions: Vec<(Box<str>, u32, u32)> = exons
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

    // Step 2. Get the reference transcript's exonic regions
    let reference_exons: Vec<(Box<str>, u32, u32)> = reference_transcript
        .get_sorted_exons()
        .iter()
        .map(|exon| (exon.chromosome.clone(), exon.start, exon.end))
        .collect();

    // Step 3. Score
    let score: f32 = match scoring_method {
        ReferenceTranscriptScoringMethod::NetOverlap => {
            let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
            let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
            let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
            num_overlap_bases - num_nonoverlap_bases
        },
        ReferenceTranscriptScoringMethod::WeightedNetOverlap => {
            let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
            let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
            let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
            num_overlap_bases - (0.5 * num_nonoverlap_bases)
        },
        ReferenceTranscriptScoringMethod::Jaccard => {
            let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
            let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
            let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
            num_overlap_bases / (num_overlap_bases + num_nonoverlap_bases)
        },
        ReferenceTranscriptScoringMethod::Overlap => {
            let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
            num_overlap_bases
        },
        ReferenceTranscriptScoringMethod::Nonoverlap => {
            let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
            let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
            num_nonoverlap_bases
        },
        ReferenceTranscriptScoringMethod::CosineSimilarity => {
            // Identify the portion of the reference gene covered by the transcript model exons
            let chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_transcript.chromosome).unwrap();
            let mut exon_reference_start_positions: Vec<u32> = Vec::new();
            let mut exon_reference_end_positions: Vec<u32> = Vec::new();
            let reference_gene_start: isize = reference_gene.start as isize;
            let reference_gene_end: isize = reference_gene.end as isize;
            for exon in exons.iter() {
                if exon.reference_chromosome_id == chromosome_id {
                    if overlaps(exon.reference_start as isize, exon.reference_end as isize, reference_gene_start, reference_gene_end) {
                        exon_reference_start_positions.push(exon.reference_start);
                        exon_reference_end_positions.push(exon.reference_end);
                    }
                }
            }
            if exon_reference_start_positions.is_empty() || exon_reference_end_positions.is_empty() {
                // Transcript model exons do not overlap the reference gene region
                0.0f32
            } else {
                // Transcript model exons overlap the reference gene region
                let min_exon_reference_start_position = *exon_reference_start_positions.iter().min().expect("Exon reference start positions vector is empty.");
                let max_exon_reference_end_position = *exon_reference_end_positions.iter().max().expect("Exon reference end positions vector is empty.");
                
                // Tight boundaries
                let vectorization_reference_start: u32 = reference_transcript.start.max(min_exon_reference_start_position);
                let vectorization_reference_end: u32 = reference_transcript.end.min(max_exon_reference_end_position);

                assert!(
                    vectorization_reference_start <= vectorization_reference_end, 
                    "vectorization_reference_start ({}) is expected to be smaller or equal to vectorization_reference_end ({})",
                    vectorization_reference_start, vectorization_reference_end
                );

                let vectorized_transcript: Vec<i8> = vectorize_exons(
                    exons,
                    chromosome_id,
                    vectorization_reference_start,
                    vectorization_reference_end,
                    1 as i8,
                    0 as i8
                );

                let vectorized_reference_transcript: Vec<i8> = reference_transcript.vectorize_exons(
                    reference_transcript.chromosome.clone(),
                    vectorization_reference_start,
                    vectorization_reference_end,
                    1 as i8,
                    0 as i8
                );

                calculate_cosine_similarity(&vectorized_transcript, &vectorized_reference_transcript) as f32
            }
        },
        ReferenceTranscriptScoringMethod::L2Distance => {
            let chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_transcript.chromosome).unwrap();
            let vectorized_transcript: Vec<i8> = vectorize_exons(
                exons,
                chromosome_id,
                reference_transcript.start,
                reference_transcript.end,
                1 as i8,
                0 as i8
            );
            let vectorized_reference_transcript: Vec<i8> = reference_transcript.vectorize_exons(
                reference_transcript.chromosome.clone(),
                reference_transcript.start,
                reference_transcript.end,
                1 as i8,
                0 as i8
            );
            calculate_l2_distance(&vectorized_transcript, &vectorized_reference_transcript) as f32
        },
        _ => {
            panic!("Unexpected scoring method: {}", scoring_method.as_str());
        }
    };

    let num_overlap_bases: u32 = count_common_bases(&model_exon_regions, &reference_exons);

    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(
        &model_exon_regions,
        &reference_exons
    );

    let reference_transcript_match: ReferenceTranscriptMatch = ReferenceTranscriptMatch::new(
        &*reference_transcript.gene_id,
        &*reference_gene.gene_name,
        &*reference_transcript.transcript_id,
        num_overlap_bases,
        num_model_only_bases,
        num_reference_only_bases,
        scoring_method.clone(),
        score
    );

    reference_transcript_match
}

fn score_reference_transcripts(
    exons: &Vec<AssembledTranscriptExon>,
    reference_transcripts: Vec<&Transcript>,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    scoring_method: ReferenceTranscriptScoringMethod
) -> Vec<ReferenceTranscriptMatch> {
    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = reference_transcripts
        .iter()
        .map(|reference_transcript| {
            let reference_transcript_match: ReferenceTranscriptMatch = score_reference_transcript(
                exons,
                reference_transcript,
                reference_gene,
                chromosome_names_map,
                scoring_method.clone()
            );
            reference_transcript_match
        })
        .collect::<Vec<_>>();
    reference_transcript_matches
}
