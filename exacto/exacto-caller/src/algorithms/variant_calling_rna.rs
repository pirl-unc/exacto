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
use exacto_core::log_info;
use indicatif::{ProgressBar, ProgressStyle};
use noodles_bam as bam;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use tempfile::{NamedTempFile, TempPath};

use crate::prelude::*;


pub fn identify_variant_transcripts(
    bam_file: &str,
    bam_bai_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    reference_transcript_scoring_method: ReferenceTranscriptScoringMethod,
    reference_transcript_selection_strategy: ReferenceTranscriptSelectionStrategy,
    top_k: usize,
    threshold: f32,
    min_mapping_quality: u16,
    min_base_quality: u8,
    num_threads: usize,
    chunk_size: usize,
    temp_dir: &str
) -> TranscriptModelSet {
    // Step 1. Get a map of read names and IDs
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        num_threads
    );

    // Step 2. Get chromosome names map
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 3. Fetch all BAM records
    log_info!("Fetching all BAM records");
    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        num_threads
    );

    // Step 4. Fetch the temp directory
    let dir: String = if temp_dir.is_empty() {
        // Use TMPDIR environment variable or fallback to system temp directory
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }
    let mut temp_files: Vec<TempPath> = Vec::new();

    // Step 5. Construct transcript models
    log_info!("Constructing transcript models");
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let entries: Vec<(usize, Vec<bam::Record>)> = records_map.into_iter().collect();
    log_info!("{} transcript models to process", entries.len());
    for chunk in entries.chunks(chunk_size) {
        let assembled_transcripts: Vec<(AssembledTranscript, Vec<TranscriptModel>)> = thread_pool.install(|| {
            chunk
                .into_par_iter()
                .filter_map(|(read_id, records)| {
                    let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();

                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        (|| -> Result<Option<(AssembledTranscript, Vec<TranscriptModel>)>, Box<dyn std::error::Error + Send + Sync>> {
                            let read_sequence: Box<str> = get_fastx_read_sequence(records);
                            let base_quality_scores: Vec<u8> = get_fastx_base_quality_scores(records);

                            let alignment: Alignment = Alignment::new(
                                *read_id,
                                &*read_sequence,
                                &base_quality_scores,
                                &records.iter().map(|record| Arc::new(record.clone())).collect()
                            );

                            let max_mapping_quality = alignment.get_alignment_records().iter()
                                .map(|r| r.record.mapping_quality().unwrap().get() as u16)
                                .max()
                                .unwrap_or(0);

                            if min_mapping_quality > max_mapping_quality {
                                return Ok(None);
                            }

                            let _ = *read_id; // read_id retained for downstream BiMap key compatibility
                            let assembled_transcript = AssembledTranscript::new(
                                &*read_name,
                                alignment.get_alignment_structure(),
                                &chromosome_names_map,
                                reference_genome_fasta_file
                            );

                            let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
                                assembled_transcript.get_exons(),
                                gene_annotator,
                                &chromosome_names_map,
                                reference_transcript_scoring_method.clone(),
                                reference_transcript_selection_strategy.clone(),
                                top_k,
                                threshold
                            );

                            let transcript_models: Vec<TranscriptModel> = assembled_transcript.identify_variants(
                                &reference_transcript_matches,
                                gene_annotator,
                                reference_genome_fasta_file,
                                &chromosome_names_map,
                                min_mapping_quality,
                                min_base_quality
                            );

                            Ok(Some((assembled_transcript, transcript_models)))
                        })()
                    }));

                    match result {
                        Ok(Ok(maybe_pair)) => maybe_pair,
                        Ok(Err(e)) => {
                            eprintln!("Error processing read name {}: {}", read_name, e);
                            None
                        }
                        Err(_panic) => {
                            eprintln!("Panic while processing read name {}", read_name);
                            None
                        }
                    }
                })
                .collect()
        });

        log_info!("\t{} assembled transcript(s) constructed in chunk.", assembled_transcripts.len());

        // Serialize chunk to temp file
        let temp_path: TempPath = {
            let temp_file = NamedTempFile::new_in(dir_path).unwrap();
            let mut writer = BufWriter::new(temp_file);
            bincode::serialize_into(&mut writer, &assembled_transcripts)
                .expect("Failed to serialize assembled transcripts");
            writer.flush().unwrap();
            let temp_file: NamedTempFile = writer.into_inner().unwrap();
            temp_file.into_temp_path()
        };
        temp_files.push(temp_path);
    }

    // Step 6. Load temp files and merge into TranscriptModelSet
    log_info!("Loading temp files and merging into a transcript model set.");
    let transcript_model_set: TranscriptModelSet = merge_temp_transcript_model_sets(
        &temp_files,
        &read_names_map,
        &chromosome_names_map,
        num_threads
    );

    transcript_model_set
}

fn merge_temp_transcript_model_sets(
    temp_files: &[TempPath],
    read_names_map: &BiMap<Box<str>, usize>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    num_threads: usize
) -> TranscriptModelSet {
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let all_pairs: Vec<Vec<(AssembledTranscript, Vec<TranscriptModel>)>> = thread_pool.install(|| {
        temp_files
            .par_iter()
            .map(|temp_path| {
                let file = File::open(temp_path).unwrap();
                let reader = BufReader::new(file);
                bincode::deserialize_from(reader)
                    .expect("Failed to deserialize assembled transcripts")
            })
            .collect()
    });

    let mut transcript_model_set = TranscriptModelSet::new();
    let mut next_transcript_model_id: usize = 1;
    for pairs in all_pairs {
        for (assembled_transcript, mut transcript_models) in pairs {
            // Only keep variant transcript models
            let keep: bool = transcript_models
                .iter()
                .any(|tm| !tm.is_reference_transcript());
            if !keep {
                continue;
            }

            for transcript_model in transcript_models.iter_mut() {
                transcript_model.set_id(next_transcript_model_id);
                next_transcript_model_id += 1;
            }

            transcript_model_set.add_transcript_model(assembled_transcript, transcript_models);
        }
    }

    transcript_model_set.load_read_names(read_names_map.clone());
    transcript_model_set.load_chromosome_names(chromosome_names_map.clone());

    transcript_model_set
}