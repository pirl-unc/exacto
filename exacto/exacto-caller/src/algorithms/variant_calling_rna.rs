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
use indicatif::{ProgressBar, ProgressStyle};
use noodles_bam as bam;
use rayon::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::log_info;
use crate::macros::*;
use crate::prelude::*;
use crate::structs::reference_transcript_match::ReferenceTranscriptMatch;
use crate::structs::transcript_model_set::TranscriptModelSet;


pub fn score_reference_transcript(
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
                .get_by_right(&exon.chromosome_id)
                .unwrap()
                .to_string()
                .into_boxed_str();
            (chr, exon.start, exon.end)
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
            let mut transcript_model: TranscriptModel = TranscriptModel::new(
                0,
                Vec::new(),
                Vec::new()
            );
            for exon in exons.iter() {
                transcript_model.add_exon(exon.clone());
            }

            // Identify the portion of the reference gene covered by the transcript model exons
            let chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_transcript.chromosome).unwrap();
            let mut start_positions: Vec<u32> = Vec::new();
            let mut end_positions: Vec<u32> = Vec::new();
            let reference_gene_start: isize = reference_gene.start as isize;
            let reference_gene_end: isize = reference_gene.end as isize;
            for exon in transcript_model.exons.iter() {
                if exon.chromosome_id == chromosome_id {
                    if overlaps(exon.start as isize, exon.end as isize, reference_gene_start, reference_gene_end) {
                        start_positions.push(exon.start);
                        end_positions.push(exon.end);
                    }
                }
            }
            if start_positions.is_empty() || end_positions.is_empty() {
                // Transcript model exons do not overlap the reference gene region
                0.0f32
            } else {
                // Transcript model exons overlap the reference gene region
                let min_start_position = *start_positions.iter().min().expect("Exon start positions vector is empty.");
                let max_end_position = *end_positions.iter().max().expect("Exon end positions vector is empty.");
                let start: u32 = reference_transcript.start.min(min_start_position);
                let end: u32 = reference_transcript.end.max(max_end_position);

                assert!(start <= end, "start ({}) is expected to be smaller or equal to end ({})", start, end);

                let vectorized_transcript: Vec<i8> = transcript_model.vectorize_exons(
                    chromosome_id,
                    start,
                    end,
                    1 as i8,
                    0 as i8
                );
                let vectorized_reference_transcript: Vec<i8> = reference_transcript.vectorize_exons(
                    reference_transcript.chromosome.clone(),
                    start,
                    end,
                    1 as i8,
                    0 as i8
                );
                calculate_cosine_similarity(&vectorized_transcript, &vectorized_reference_transcript) as f32
            }
        },
        ReferenceTranscriptScoringMethod::L2Distance => {
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
            calculate_l2_distance(&vectorized_transcript, &vectorized_reference_transcript) as f32
        },
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
        scoring_method.clone(),
        score
    );

    reference_transcript_match
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

pub fn identify_reference_transcript_matches(
    exons: &Vec<TranscriptModelExon>,
    reference_genes: Vec<&Gene>,
    chromosome_names_map: &BiMap<Box<str>,u16>,
    scoring_method: ReferenceTranscriptScoringMethod,
    selection_strategy: ReferenceTranscriptSelectionStrategy,
    top_k: usize,
    threshold: f32
) -> Vec<ReferenceTranscriptMatch> {
    // Step 1. Score each transcript in each gene
    let mut reference_transcript_scores: HashMap<Box<str>,Vec<ReferenceTranscriptMatch>> = HashMap::new();
    for reference_gene in reference_genes.iter() {
        let mut reference_transcripts: Vec<&Transcript> = Vec::new();
        for reference_transcript in reference_gene.transcripts.values() {
            reference_transcripts.push(reference_transcript);
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
    
    // Step 2. Select reference transcript matches
    let mut reference_transcript_matches: Vec<ReferenceTranscriptMatch> = Vec::new();
    for reference_gene in reference_genes.iter() {
        if let Some(matches) = reference_transcript_scores.get(&reference_gene.gene_id) {
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

    reference_transcript_matches
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
    reference_transcript_scoring_method: ReferenceTranscriptScoringMethod,
    reference_transcript_selection_strategy: ReferenceTranscriptSelectionStrategy,
    top_k: usize,
    threshold: f32,
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
    log_info!("Fetching all BAM records");
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        num_threads
    );

    // Step 4. Construct transcript models
    log_info!("Constructing transcript models");
    let pb = Arc::new(ProgressBar::new(records_map.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let transcript_models: Vec<TranscriptModel> = thread_pool.install(|| {
        records_map
            .par_iter()
            .filter_map(|(read_id, records)| {
                let result = {
                    // Get read original sequence
                    let read_sequence: Box<str> = get_original_read_sequence(records.iter().collect::<Vec<_>>().as_slice());

                    // Get base quality scores
                    let quality_scores: Vec<u8> = get_original_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());

                    // Construct an Alignment object
                    let alignment: Alignment = Alignment::new(
                        *read_id,
                        read_sequence,
                        quality_scores,
                        records.clone(),
                    );

                    // Get the maximum mapping quality score
                    let mut min_mapping_quality_ = 0u8;
                    for alignment_record in alignment.alignment_records.iter() {
                        let curr_mapping_quality = alignment_record.record.mapping_quality().unwrap().get();
                        min_mapping_quality_ = min_mapping_quality_.max(curr_mapping_quality);
                    }
                    if min_mapping_quality > min_mapping_quality_ as usize {
                        return None;
                    }

                    // Identify exons
                    let exons = alignment.identify_exons(min_mapping_quality);

                    // Identify splice junctions
                    let splice_junctions = alignment.identify_splice_junctions(
                        &chromosome_names_map,
                        reference_genome_fasta_file,
                        min_mapping_quality,
                    );

                    // Identify reference gene IDs
                    let reference_gene_ids = identify_overlapping_gene_ids(
                        &exons,
                        gene_annotator,
                        &chromosome_names_map,
                    );

                    // Fetch reference genes
                    let mut reference_genes: Vec<&Gene> = Vec::new();
                    for reference_gene_id in reference_gene_ids.iter() {
                        let reference_gene = gene_annotator.get_gene(reference_gene_id).unwrap();
                        reference_genes.push(reference_gene);
                    }

                    // Identify closest reference transcripts
                    let reference_transcript_matches = identify_reference_transcript_matches(
                        &exons,
                        reference_genes,
                        &chromosome_names_map,
                        reference_transcript_scoring_method.clone(),
                        reference_transcript_selection_strategy.clone(),
                        top_k,
                        threshold
                    );

                    // Identify splice variants
                    let splice_variant_records: HashMap<Box<str>,Vec<VariantRecord>> = alignment.identify_splice_variant_records(
                        &reference_transcript_matches,
                        &chromosome_names_map,
                        min_mapping_quality,
                    );

                    // Identify sequence variants
                    let sequence_variant_records = alignment.identify_sequence_variant_records(
                        min_mapping_quality,
                        min_average_base_quality,
                    );

                    // Construct a TranscriptModel object
                    let mut transcript_model = TranscriptModel::new(
                        alignment.read_id,
                        reference_transcript_matches,
                        vec![*read_id],
                    );
                    for exon in exons {
                        transcript_model.add_exon(exon);
                    }
                    for splice_junction in splice_junctions {
                        transcript_model.add_splice_junction(splice_junction);
                    }
                    for (reference_transcript_id, variant_records) in splice_variant_records.iter() {
                        for variant_record in variant_records.iter() {
                            let mut variant_call = VariantCall::new();
                            variant_call.add_variant_record(variant_record.clone());
                            transcript_model.add_splice_variant_call(&*reference_transcript_id, variant_call);
                        }
                    }
                    for variant_record in sequence_variant_records {
                        let mut variant_call = VariantCall::new();
                        variant_call.add_variant_record(variant_record);
                        transcript_model.add_sequence_variant_call(variant_call);
                    }

                    Some(transcript_model)
                };
                pb.inc(1);
                result
            })
            .collect()
    });
    pb.finish_with_message("Completed constructing transcript models.");

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
