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


fn identify_overlapping_gene_ids(
    exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>
) -> HashSet<Box<str>> {
    let mut gene_ids: HashSet<Box<str>> = HashSet::new();
    for exon in exons.iter() {
        let chromosome_name: Box<str> = chromosome_names_map.get_by_right(&exon.reference_chromosome_id).unwrap().to_string().into_boxed_str();
        let gene_ids_: Vec<Box<str>> = gene_annotator.get_gene_ids_overlapping_region(&*chromosome_name, exon.reference_start, exon.reference_end);
        for gene_id in gene_ids_ {
            let gene: &Gene = gene_annotator.get_gene(&*gene_id).unwrap();
            if gene.strand.as_str() == exon.reference_strand.as_str() {
                gene_ids.insert(gene_id);
            }
        }
    }
    gene_ids
}

fn identify_overlapping_transcript_ids(
    exons: &Vec<TranscriptModelExon>,
    reference_gene_ids: &HashSet<Box<str>>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>
) -> HashSet<Box<str>> {
    let mut transcript_ids: HashSet<Box<str>> = HashSet::new();
    for reference_gene_id in reference_gene_ids {
        let reference_gene: &Gene = gene_annotator.get_gene(&reference_gene_id).unwrap();
        for transcript in reference_gene.transcripts.values() {
            let transcript: &Transcript = gene_annotator.get_transcript(&*transcript.transcript_id).unwrap();
            for exon in exons.iter() {
                let chromosome_name: Box<str> = chromosome_names_map.get_by_right(&exon.reference_chromosome_id).unwrap().to_string().into_boxed_str();
                if transcript.chromosome == chromosome_name &&
                    transcript.strand == exon.reference_strand {
                    if overlaps(
                        exon.reference_start as isize,
                        exon.reference_end as isize,
                        transcript.start as isize,
                        transcript.end as isize
                    ) {
                        transcript_ids.insert(transcript.transcript_id.clone());
                    }
                }
            }
        }
    }
    transcript_ids
}

fn score_reference_transcript(
    exons: &Vec<TranscriptModelExon>,
    reference_transcript: &Transcript,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<Box<str>,u16>,
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

                // Loose boundaries 
                // let vectorization_reference_start: u32 = reference_transcript.start.min(min_exon_reference_start_position);
                // let vectorization_reference_end: u32 = reference_transcript.end.max(max_exon_reference_end_position);

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

    let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;

    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);

    let reference_transcript_match: ReferenceTranscriptMatch = ReferenceTranscriptMatch::new(
        &*reference_transcript.gene_id,
        &*reference_transcript.transcript_id,
        num_overlap_bases as u32,
        num_model_only_bases,
        num_reference_only_bases,
        scoring_method.clone(),
        score
    );

    reference_transcript_match
}

pub fn identify_reference_transcript_matches(
    exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>,
    scoring_method: ReferenceTranscriptScoringMethod,
    selection_strategy: ReferenceTranscriptSelectionStrategy,
    top_k: usize,
    threshold: f32
) -> Vec<ReferenceTranscriptMatch> {
    // Step 1. Identifying overlapping genes and transcript IDs
    let reference_gene_ids: HashSet<Box<str>> = identify_overlapping_gene_ids(
        exons,
        gene_annotator,
        chromosome_names_map
    );
    let reference_transcript_ids: HashSet<Box<str>> = identify_overlapping_transcript_ids(
        exons,
        &reference_gene_ids,
        gene_annotator,
        chromosome_names_map
    );

    // Step 2. Score each transcript in each gene
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

    // Step 3. Select reference transcript matches
    let mut reference_transcript_matches: Vec<ReferenceTranscriptMatch> = Vec::new();
    for reference_gene_id in reference_gene_ids.iter() {
        if let Some(matches) = reference_transcript_scores.get(reference_gene_id) {
            // Sort by score descending
            let mut sorted_matches = matches.clone();
            sorted_matches.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

            if selection_strategy == ReferenceTranscriptSelectionStrategy::TopK {
                // Find the top K score values
                let mut prev_score: f32 = f32::NEG_INFINITY;
                let mut top_k_scores: Vec<f32> = Vec::new();
                for reference_transcript_match in sorted_matches.iter() {
                    let score = reference_transcript_match.score;
                    if (score == prev_score) || (top_k_scores.len() >= top_k) {
                        continue;
                    }
                    top_k_scores.push(score);
                    prev_score = score;
                }

                // Find reference transcript matches that have the top K score values
                let min_score = top_k_scores
                    .iter()
                    .copied()
                    .reduce(|a, b| a.partial_cmp(&b).map(|o| if o == std::cmp::Ordering::Less { a } else { b }).unwrap())
                    .expect("top_k_scores is empty");
                for reference_transcript_match in sorted_matches.iter() {
                    if reference_transcript_match.score >= min_score {
                        reference_transcript_matches.push(reference_transcript_match.clone());
                    }
                }
            } else if selection_strategy == ReferenceTranscriptSelectionStrategy::Threshold {
                reference_transcript_matches.extend(
                    sorted_matches.into_iter().filter(|m| m.score >= threshold)
                );
            } else {
                panic!("Unsupported selection strategy: {}", selection_strategy.as_str());
            }
        }
    }

    // Step 4. Sort by score
    reference_transcript_matches.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    reference_transcript_matches
}

pub fn score_reference_transcripts(
    exons: &Vec<TranscriptModelExon>,
    reference_transcripts: Vec<&Transcript>,
    reference_gene: &Gene,
    chromosome_names_map: &BiMap<Box<str>,u16>,
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
