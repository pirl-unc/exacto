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
use exacto_util::prelude::*;
use exacto_util as util;
use noodles_bam as bam;
use rayon::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use std::collections::{HashMap, HashSet};

use crate::prelude::*;
use crate::structs::reference_transcript_match::ReferenceTranscriptMatch;
use crate::structs::transcript_model_set::TranscriptModelSet;


pub fn identify_closest_reference_transcript_ids_helper(
    exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>,
    gene_id: &str,
    scoring_method: ReferenceTranscriptScoringMethods
) -> ReferenceTranscriptMatch {
    // Step 1. Get the transcript model's exonic regions
    let model_exon_regions: Vec<(Box<str>, u32, u32)> = exons
        .iter()
        .map(|exon| {
            let chr = chromosome_names_map
                .get_by_right(&exon.chromosome_id)
                .unwrap()
                .to_string()
                .into_boxed_str();
            (chr, exon.start, exon.end)
        })
        .collect();

    // Step 2. Score each reference transcript
    let reference_gene = gene_annotator.get_gene(gene_id).unwrap();
    let mut matched_reference_transcripts: Vec<ReferenceTranscriptMatch> = reference_gene
        .get_transcript_ids()
        .iter()
        .map(|reference_transcript_id| {
            let reference_transcript = gene_annotator.get_transcript(&**reference_transcript_id).unwrap();
            let reference_exons: Vec<(Box<str>, u32, u32)> = reference_transcript
                .get_sorted_exons()
                .iter()
                .map(|exon| (exon.chromosome.clone(), exon.start, exon.end))
                .collect();

            let score: f32 = match scoring_method {
                ReferenceTranscriptScoringMethods::NetOverlap => {
                    let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
                    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
                    let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
                    num_overlap_bases - num_nonoverlap_bases
                },
                ReferenceTranscriptScoringMethods::WeightedNetOverlap => {
                    let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
                    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
                    let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
                    num_overlap_bases - (0.5 * num_nonoverlap_bases)
                },
                ReferenceTranscriptScoringMethods::Jaccard => {
                    let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
                    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
                    let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
                    num_overlap_bases / (num_overlap_bases + num_nonoverlap_bases)
                },
                ReferenceTranscriptScoringMethods::Overlap => {
                    let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
                    num_overlap_bases
                },
                ReferenceTranscriptScoringMethods::Nonoverlap => {
                    let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);
                    let num_nonoverlap_bases: f32 = num_model_only_bases as f32 + num_reference_only_bases as f32;
                    num_nonoverlap_bases
                },
                ReferenceTranscriptScoringMethods::CosineSimilarity => {
                    let mut transcript_model: TranscriptModel = TranscriptModel::new(
                        0,
                        Vec::new(),
                        Vec::new()
                    );
                    for exon in exons.iter() {
                        transcript_model.add_exon(exon.clone());
                    }
                    let chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_transcript.chromosome).unwrap();
                    let vectorized_transcript: Vec<i8> = transcript_model.vectorize_exons(
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
                    calculate_cosine_similarity(&vectorized_transcript, &vectorized_reference_transcript) as f32
                }
                _ => {
                    panic!("Unexpected scoring method: {}", scoring_method.as_str());
                }
            };

            let num_overlap_bases: f32 = count_common_bases(&model_exon_regions, &reference_exons) as f32;
            let (num_model_only_bases, num_reference_only_bases) = count_non_overlapping_bases(&model_exon_regions, &reference_exons);

            let reference_transcript_match: ReferenceTranscriptMatch = ReferenceTranscriptMatch::new(
                reference_transcript.clone(),
                num_overlap_bases as u32,
                num_model_only_bases,
                num_reference_only_bases,
                score
            );

            reference_transcript_match
        })
        .collect::<Vec<_>>();

    // Step 3. Identify closest transcript
    if matched_reference_transcripts.is_empty() {
        let transcript_: Transcript = Transcript::new(
            "",
            "",
            "",
            "",
            0,
            0,
            util::common::constants::Strands::Forward,
            0,
            "",
            "",
            ""
        );
        let reference_transcript_match: ReferenceTranscriptMatch = ReferenceTranscriptMatch::new(
            transcript_,
            0u32,
            0u32,
            0u32,
            0.0f32
        );
        return reference_transcript_match;
    }

    // Get the max score
    let best_match = matched_reference_transcripts
        .into_iter()
        .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap());

    best_match.unwrap()
}

pub fn identify_closest_reference_transcript_ids(
    exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>,
    gene_ids: &HashSet<Box<str>>,
    scoring_method: ReferenceTranscriptScoringMethods
) -> Vec<ReferenceTranscriptMatch> {
    // Identify the closest reference transcript ID for each gene
    let mut reference_transcript_scores: HashMap<Box<str>,ReferenceTranscriptMatch> = HashMap::new();
    for reference_gene_id in gene_ids.iter() {
        let matched_reference_transcript: ReferenceTranscriptMatch = identify_closest_reference_transcript_ids_helper(
            exons,
            gene_annotator,
            chromosome_names_map,
            &*reference_gene_id,
            scoring_method.clone()
        );
        reference_transcript_scores.insert(reference_gene_id.clone(),matched_reference_transcript);
    }

    // Choose the closest reference transcript ID if genes overlap
    let mut reference_gene_ids: HashSet<Box<str>> = HashSet::new();
    for reference_gene_id_1 in reference_transcript_scores.keys() {
        let reference_gene_1 = gene_annotator.get_gene(reference_gene_id_1).unwrap();
        let (start_1, end_1) = (reference_gene_1.start as isize, reference_gene_1.end as isize);

        // Collect all overlapping gene IDs, including itself
        let overlapping_gene_ids: HashSet<_> = reference_transcript_scores
            .keys()
            .filter(|reference_gene_id_2| {
                if *reference_gene_id_1 == **reference_gene_id_2 {
                    true
                } else {
                    let reference_gene_2 = gene_annotator.get_gene(reference_gene_id_2).unwrap();
                    overlaps(start_1, end_1, reference_gene_2.start as isize, reference_gene_2.end as isize)
                }
            })
            .cloned()
            .collect();

        // Pick the best-scoring gene from the overlapping set
        if let Some(best_gene_id) = overlapping_gene_ids
            .iter()
            .max_by(|a, b| {
                let matched_reference_transcript_a: &ReferenceTranscriptMatch = reference_transcript_scores.get(*a).unwrap();
                let matched_reference_transcript_b: &ReferenceTranscriptMatch = reference_transcript_scores.get(*b).unwrap();
                let score_a = matched_reference_transcript_a.score;
                let score_b = matched_reference_transcript_b.score;
                score_a.partial_cmp(&score_b).unwrap()
            }) {
            reference_gene_ids.insert(best_gene_id.clone());
        }
    }

    // Retain the best scoring gene transcripts
    let reference_transcript_ids: Vec<ReferenceTranscriptMatch> = reference_gene_ids
        .iter()
        .filter_map(|gene_id| reference_transcript_scores.get(gene_id).map(|reference_transcript_match| reference_transcript_match.clone()))
        .collect();

    reference_transcript_ids
}

pub fn identify_overlapping_gene_ids(
    exons: &Vec<TranscriptModelExon>,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>
) -> HashSet<Box<str>> {
    let mut gene_ids: HashSet<Box<str>> = HashSet::new();
    for exon in exons.iter() {
        let chromosome: Box<str> = chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap().to_string().into_boxed_str();
        let gene_ids_: Vec<Box<str>> = gene_annotator.get_gene_ids_overlapping_region(&*chromosome, exon.start, exon.end);
        for gene_id in gene_ids_ {
            let gene: &Gene = gene_annotator.get_gene(&*gene_id).unwrap();
            if gene.strand.as_str() == exon.strand.as_str() {
                gene_ids.insert(gene_id);
            }
        }
    }
    gene_ids
}

pub fn identify_variant_transcripts(
    bam_file: &str,
    bam_bai_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    reference_transcript_scoring_method: ReferenceTranscriptScoringMethods,
    min_mapping_quality: usize,
    min_average_base_quality: f32,
    num_threads: usize
) -> TranscriptModelSet {
    // Step 1. Get a map of read names and IDs
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        num_threads
    );

    // Step 2. Get chromosome names map
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);

    // Step 3. Fetch all BAM records
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        num_threads
    );

    // Step 4. Construct transcript models
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let transcript_models: Vec<TranscriptModel> = thread_pool.install(|| {
        records_map
            .par_iter()
            .filter_map(|(read_id, records)| {
                let read_sequence: Box<str> =  get_original_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
                let quality_scores: Vec<u8> = get_original_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
                let alignment: Alignment = Alignment::new(
                    *read_id,
                    read_sequence,
                    quality_scores,
                    records.clone()
                );

                // Check mapping quality
                let mut min_mapping_quality_: u8 = 0u8;
                for alignment_record in alignment.alignment_records.iter() {
                    let curr_mapping_quality = alignment_record.record.mapping_quality().unwrap().get();
                    if min_mapping_quality_ < curr_mapping_quality {
                        min_mapping_quality_ = curr_mapping_quality;
                    }
                }
                if min_mapping_quality > min_mapping_quality_ as usize {
                    // No exon will be identified
                    return None;
                }

                // Identify exons
                let exons: Vec<TranscriptModelExon> = alignment.identify_exons(min_mapping_quality);

                // Identify splice junctions
                let splice_junctions: Vec<TranscriptModelSpliceJunction> = alignment.identify_splice_junctions(
                    &chromosome_names_map,
                    reference_genome_fasta_file,
                    min_mapping_quality
                );

                // Identify overlapping gene IDs
                let reference_gene_ids: HashSet<Box<str>> = identify_overlapping_gene_ids(
                    &exons,
                    gene_annotator,
                    &chromosome_names_map
                );

                // Identify closest reference transcript IDs
                let matched_reference_transcripts: Vec<ReferenceTranscriptMatch> = identify_closest_reference_transcript_ids(
                    &exons,
                    gene_annotator,
                    &chromosome_names_map,
                    &reference_gene_ids,
                    reference_transcript_scoring_method.clone()
                );

                // Identify splice variant records
                let splice_variant_records: Vec<VariantRecord> = alignment.identify_splice_variant_records(
                    &matched_reference_transcripts,
                    &chromosome_names_map,
                    min_mapping_quality
                );

                // Identify sequence variant records
                let sequence_variant_records: Vec<VariantRecord> = alignment.identify_sequence_variant_records(
                    min_mapping_quality,
                    min_average_base_quality
                );

                // Construct a transcript model
                let read_ids: Vec<usize> = vec![*read_id];
                let mut transcript_model: TranscriptModel = TranscriptModel::new(
                    alignment.read_id,
                    matched_reference_transcripts,
                    read_ids
                );
                for exon in exons {
                    transcript_model.add_exon(exon);
                }
                for splice_junction in splice_junctions {
                    transcript_model.add_splice_junction(splice_junction);
                }
                for splice_variant_record in splice_variant_records {
                    let mut variant_call: VariantCall = VariantCall::new();
                    variant_call.add_variant_record(splice_variant_record);
                    transcript_model.add_variant_call(variant_call);
                }
                for sequence_variant_record in sequence_variant_records {
                    let mut variant_call: VariantCall = VariantCall::new();
                    variant_call.add_variant_record(sequence_variant_record);
                    transcript_model.add_variant_call(variant_call);
                }

                Some(transcript_model)
            })
            .collect()
    });

    let mut transcript_model_set: TranscriptModelSet = TranscriptModelSet::new();
    let mut transcript_id: usize = 1;;
    for mut transcript_model in transcript_models {
        transcript_model.transcript_id = transcript_id;
        transcript_model_set.add_transcript_model(transcript_model);
        transcript_id += 1;
    }

    transcript_model_set.load_read_names(read_names_map);
    transcript_model_set.load_chromosome_names(chromosome_names_map);

    transcript_model_set
}

// pub fn is_nascent_transcript(
//     splice_junctions: &Vec<TranscriptModelSpliceJunction>,
//     exons: &Vec<TranscriptModelExon>,
//     reference_transcript_ids: Vec<&str>,
//     gene_annotator: &impl GeneAnnotator
// ) -> bool {
//     if splice_junctions.is_empty() {
//         // The transcript model is a nascent transcript if it does not have any evidence of
//         // splicing and does not overlap with any known gene
//         if reference_transcript_ids.is_empty() {
//             return true;
//         }
//
//         // The transcript model is a nascent transcript if it does not have any evidence of
//         // splicing and overlaps with an intron
//         let mut overlaps_intron: bool = false;
//         for reference_transcript_id in reference_transcript_ids.iter() {
//             let reference_transcript: &Transcript = gene_annotator.get_transcript(reference_transcript_id).unwrap();
//             if reference_transcript.get_exon_ids().len() == 1 {
//                 let reference_exon: &Exon = reference_transcript.get_exon(reference_transcript.get_exon_ids().first().unwrap()).unwrap();
//                 for exon in exons.iter() {
//                     if overlaps(exon.start as isize, exon.end as isize, reference_exon.start as isize, reference_exon.end as isize) {
//                         // The transcript model is not a nascent transcript if it overlaps
//                         // with a transcript with only 1 exon
//                         return false;
//                     }
//                 }
//             } else {
//                 let mut overlaps_intron: bool = false;
//                 for reference_intron in reference_transcript.get_introns().iter() {
//                     for exon in exons.iter() {
//                         if overlaps(exon.start as isize, exon.end as isize, reference_intron.start as isize, reference_intron.end as isize) {
//                             overlaps_intron = true;
//                             break;
//                         }
//                     }
//                     if overlaps_intron {
//                         break;
//                     }
//                 }
//                 let mut overlaps_exon: bool = false;
//                 for reference_exon in reference_transcript.exons.values().iter() {
//                     for exon in exons.iter() {
//                         if overlaps(exon.start as isize, exon.end as isize, reference_exon.start as isize, reference_exon.end as isize) {
//                             overlaps_exon = true;
//                             break;
//                         }
//                     }
//                     if overlaps_exon {
//                         break;
//                     }
//                 }
//                 if overlaps_exon && overlaps_intron {
//                     // The transcript model is a nascent transcript if it does not have any
//                     // evidence of splicing and overlaps with an exon and an intron of a
//                     // transcript (with more than 1 exon).
//                     return true;
//                 }
//             }
//         }
//
//         // todo change return type to enum (classes: nascent, mature, unknown) - call it transcript stage
//         // todo consider the number of exons in a transcript
//         // todo consider other ways in which there could be no splicing and
//
//         false
//     } else {
//         // The transcript model is not a nascent transcript because there is evidence of splicing
//         false
//     }
// }
