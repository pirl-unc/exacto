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
use indicatif::{ProgressBar, ProgressStyle};
use noodles_bam as bam;
use rayon::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::prelude::*;
use crate::log_info;
use crate::macros::*;


// pub fn map_variant_records_to_reference_transcripts(
//     variant_records: Vec<&VariantRecord>,
//     reference_transcript_matches: Vec<&ReferenceTranscriptMatch>,
//     chromosome_names_map: &BiMap<Box<str>,u16>,
//     gene_annotator: &(impl GeneAnnotator + Sync)
// ) -> HashMap<Box<str>, Vec<VariantRecord>> {
//     let mut variant_records_map: HashMap<Box<str>, Vec<VariantRecord>> = HashMap::new();
//     for variant_record in variant_records {
//         let chromosome_1 = chromosome_names_map.get_by_right(&variant_record.get_chromosome_1()).unwrap().clone();
//         let chromosome_2 = chromosome_names_map.get_by_right(&variant_record.get_chromosome_2()).unwrap().clone();
//         let position_1: isize = variant_record.get_position_1() as isize;
//         let position_2: isize = variant_record.get_position_2() as isize;
//         let mut matched: bool = false;
//         for reference_transcript_match in reference_transcript_matches.iter() {
//             let reference_transcript: &Transcript = gene_annotator.get_transcript(&reference_transcript_match.reference_transcript_id).unwrap();
//             let reference_transcript_chromosome: Box<str> = reference_transcript.chromosome.clone();
//             let reference_transcript_start: isize = reference_transcript.start as isize;
//             let reference_transcript_end: isize = reference_transcript.end as isize;
//             if (reference_transcript_chromosome == chromosome_1 && overlaps(position_1, position_1, reference_transcript_start, reference_transcript_end)) ||
//                 (reference_transcript_chromosome == chromosome_2 && overlaps(position_2, position_2, reference_transcript_start, reference_transcript_end)) {
//                 variant_records_map
//                     .entry(reference_transcript_match.reference_transcript_id.clone())
//                     .or_insert_with(Vec::new)
//                     .push(variant_record.clone());
//                 matched = true;
//             }
//         }
//         if matched == false {
//             variant_records_map
//                 .entry(GenicRegion::Intergenic.as_str().to_string().into())
//                 .or_insert_with(Vec::new)
//                 .push(variant_record.clone());
//         }
//     }
//
//     variant_records_map
// }

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
    min_base_quality: u8,
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
                    let read_sequence: Box<str> = get_fastx_read_sequence(records.iter().collect::<Vec<_>>().as_slice());

                    // Get base quality scores
                    let base_quality_scores: Vec<u8> = get_fastx_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());

                    // Construct an Alignment object
                    let alignment: Alignment = Alignment::new(
                        *read_id,
                        &*read_sequence,
                        &base_quality_scores,
                        records,
                    );

                    // Check if the maximum mapping quality score is above the minimum mapping quality
                    let mut max_mapping_quality = 0usize;
                    for alignment_record in alignment.get_alignment_records().iter() {
                        let mapping_quality: usize = alignment_record.record.mapping_quality().unwrap().get() as usize;
                        max_mapping_quality = max_mapping_quality.max(mapping_quality);
                    }
                    if min_mapping_quality > max_mapping_quality {
                        return None;
                    }

                    // Construct a TranscriptModel object
                    let mut transcript_model: TranscriptModel = TranscriptModel::new(
                        1,
                        alignment.get_alignment_structure(),
                        &chromosome_names_map,
                        reference_genome_fasta_file
                    );

                    // Identify reference transcript matches
                    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
                        transcript_model.get_exons(),
                        gene_annotator,
                        &chromosome_names_map,
                        reference_transcript_scoring_method.clone(),
                        reference_transcript_selection_strategy.clone(),
                        top_k,
                        threshold
                    );

                    // Identify variants
                    transcript_model.identify_variants(
                        &reference_transcript_matches,
                        gene_annotator,
                        reference_genome_fasta_file,
                        min_mapping_quality,
                        min_base_quality
                    );

                    Some(transcript_model)
                };
                pb.inc(1);
                result
            })
            .collect()
    });
    pb.finish_with_message("Completed constructing transcript models.");

    let mut transcript_model_set: TranscriptModelSet = TranscriptModelSet::new();
    let mut transcript_id: usize = 1;
    for mut transcript_model in transcript_models {
        // Check if the transcript model is a reference transcript)
        if transcript_model.is_reference_transcript() == false {
            transcript_model.set_transcript_id(transcript_id);
            transcript_model_set.add_transcript_model(transcript_model);
            transcript_id += 1;
        }
    }
    transcript_model_set.load_read_names(read_names_map);
    transcript_model_set.load_chromosome_names(chromosome_names_map);

    transcript_model_set
}
