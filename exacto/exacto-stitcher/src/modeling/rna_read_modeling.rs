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
use exacto_caller::prelude::*;
use exacto_core::log_info;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bgzf::VirtualPosition;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::ThreadPool;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};


pub(crate) fn model_rna_read(
    read_id: usize,
    records: Vec<bam::Record>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    min_mapping_quality: u16
) -> Option<TranscriptModel> {
    // Step 1. Filter based on mapping quality.
    let max_mapping_quality: u16 = records
        .iter()
        .map(get_alignment_mapping_quality)
        .max()
        .unwrap_or(0);
    if min_mapping_quality > max_mapping_quality {
        return None;
    }

    // Step 2. Get the FASTX read sequence.
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);

    // Step 3. Get the base quality scores, defaulting to Q60 if absent.
    let base_quality_scores: Vec<u8> = {
        let scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        if scores.is_empty() {
            vec![60u8; read_sequence.len()]
        } else {
            scores
        }
    };

    // Step 4. Get references to the BAM records.
    let bam_records: Vec<Arc<bam::Record>> = records
        .into_iter()
        .map(Arc::new)
        .collect();

    // Step 5. Construct an instance of AlignmentModel.
    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &base_quality_scores,
        &bam_records
    );

    // Step 6. Identify reference transcript matches.
    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        gene_annotator,
        chromosome_names_map
    );

    // Step 7. Construct an instance of TranscriptModel.
    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        gene_annotator,
        chromosome_names_map,
        fasta_map
    );

    Some(transcript_model)
}


pub(crate) fn map_rna_reads<T: Send>(
    bam_file: &str,
    read_names_map: &BiMap<Box<str>, usize>,
    record_positions_map: &HashMap<usize, Vec<VirtualPosition>>,
    num_threads: usize,
    f: impl Fn(usize, Vec<bam::Record>) -> T + Sync
) -> Vec<T> {
    let num_reads: usize = record_positions_map.len();
    let num_characterized: AtomicUsize = AtomicUsize::new(0);
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut read_ids: Vec<usize> = read_names_map.right_values().copied().collect();
    read_ids.sort();
    thread_pool.install(|| {
        read_ids
            .par_iter()
            .map_init(
                // One reader per worker thread - reused across every read ID it handles.
                || {
                    let mut reader = bam::io::reader::Builder::default()
                        .build_from_path(bam_file)
                        .unwrap();
                    reader.read_header().unwrap(); // consume header before seeking
                    reader
                },
                |reader, &read_id| -> T {
                    // Get the BAM records.
                    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(
                        reader,
                        read_id,
                        &record_positions_map
                    );

                    let result: T = f(read_id, records);

                    // Log the progress.
                    let num_done: usize = num_characterized.fetch_add(1, Ordering::Relaxed) + 1;
                    if (num_done * 10) / num_reads != ((num_done - 1) * 10) / num_reads {
                        log_info!(
                                "Characterized {}% of reads ({}/{})",
                                (num_done * 100) / num_reads,
                                num_done,
                                num_reads
                            );
                    }

                    result
                }
            )
            .collect()
    })
}


#[cfg(test)]
pub(crate) fn model_rna_reads(
    bam_file: &str,
    read_names_map: &BiMap<Box<str>, usize>,
    record_positions_map: &HashMap<usize, Vec<VirtualPosition>>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize,
    min_mapping_quality: u16
) -> Vec<TranscriptModel> {
    map_rna_reads(
        bam_file,
        read_names_map,
        record_positions_map,
        num_threads,
        |read_id, records| std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| model_rna_read(
            read_id,
            records,
            gene_annotator,
            chromosome_names_map,
            fasta_map,
            min_mapping_quality
        ))).ok().flatten()
    )
        .into_iter()
        .flatten()
        .collect()
}
