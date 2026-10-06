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
use exacto_cluster::prelude::*;
use exacto_core::log_info;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::Writer;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::prelude::*;
use crate::correction::corrected_rna_read::CorrectedRNARead;
use crate::correction::rna_read_correction::correct_rna_read;


pub fn correct_rna_reads(
    bam_file: &str,
    output_fastq_file: &str,
    cluster_set: &RNAReadClusterSet,
    gene_annotator: Option<&(dyn GeneAnnotator + Sync)>,
    options: &CorrectRNAReadsOptions,
    num_threads: usize,
    chunk_size: usize
) -> Result<(), CorrectionError> {
    let gene_annotator: Option<&(dyn GeneAnnotator + Sync)> = match (options.trim_transcript_ends, gene_annotator) {
        (false, _) => None,
        (true, None) => return Err(CorrectionError::NoGeneAnnotation),
        (true, gene_annotator) => gene_annotator
    };

    // Step 1. Index the BAM file.
    log_info!("Indexing the BAM file.");
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(
        bam_file,
        true,
        num_threads
    );

    // Step 2. Create a BiMap of chromosome names and IDs, in BAM header order.
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 3. Get cluster IDs and read names, in the order the output is written. A read is
    // corrected once, in the cluster it was assigned to: a cluster it is shared with fits it as
    // well, and a second copy in the FASTQ would be a second read.
    // Vec<(cluster ID, read name)>
    let mut read_names: Vec<(usize, Box<str>)> = cluster_set
        .get_read_clusters_map()
        .into_iter()
        .map(|(read_name, cluster_id)| (cluster_id, read_name))
        .collect();
    read_names.sort_unstable();
    let num_reads: usize = read_names.len();

    // Every read must be in the BAM before any is corrected: a wrong or truncated BAM is an
    // input error, not something to find out from a worker halfway through the output.
    let missing: Vec<&Box<str>> = read_names
        .iter()
        .map(|(_, read_name)| read_name)
        .filter(|read_name| !read_names_map.contains_left(*read_name))
        .collect();
    if let Some(&read_name) = missing.first() {
        return Err(CorrectionError::ReadsNotInBam {
            bam_file: bam_file.into(),
            num_missing: missing.len(),
            num_reads: num_reads,
            read_name: read_name.clone()
        });
    }

    // Step 4. Get each cluster's called operations, keyed by the BAM's chromosome IDs.
    let cluster_operations: HashMap<usize, HashSet<GraphOperation>> = translate_cluster_operations(
        cluster_set,
        &chromosome_names_map
    );
    let cluster_operation_refs: HashMap<usize, HashSet<&GraphOperation>> = cluster_operations
        .iter()
        .map(|(&cluster_id, operations)| (cluster_id, operations.iter().collect()))
        .collect();

    // Step 5. Cut the reads into chunks.
    let mut chunks: Vec<&[(usize, Box<str>)]> = Vec::new();
    let mut remaining: &[(usize, Box<str>)] = &read_names;
    while !remaining.is_empty() {
        let last_cluster: usize = remaining[chunk_size.clamp(1, remaining.len()) - 1].0;
        let end: usize = remaining.partition_point(|&(cluster_id, _)| cluster_id <= last_cluster);
        let (chunk, rest) = remaining.split_at(end);
        chunks.push(chunk);
        remaining = rest;
    }

    // Step 6. Correct one chunk at a time, writing as we go.
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let partial_fastq_file: String = format!("{output_fastq_file}.partial");
    let mut writer: Writer<File> = open_fastq_writer(&partial_fastq_file);
    let num_corrected: AtomicUsize = AtomicUsize::new(0);
    let num_uncorrected: AtomicUsize = AtomicUsize::new(0);
    for chunk in chunks.iter() {
        let corrected_reads: Vec<CorrectedRNARead> = thread_pool.install(|| {
            chunk
                .par_iter()
                .map_init(
                    || {
                        let mut reader = bam::io::reader::Builder::default()
                            .build_from_path(bam_file)
                            .unwrap();
                        reader.read_header().unwrap();
                        reader
                    },
                    |reader, (cluster_id, read_name)| {
                        let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
                        let records: Vec<bam::Record> = fetch_bam_records_for_read_id(reader, read_id, &record_positions_map);
                        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
                        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
                        let curr_cluster_graph_operations: &HashSet<&GraphOperation> = cluster_operation_refs
                            .get(cluster_id)
                            .unwrap();

                        let corrected_rna_read: CorrectedRNARead = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            let records: Vec<Arc<bam::Record>> = records
                                .into_iter()
                                .map(Arc::new)
                                .collect();
                            let alignment_model: AlignmentModel = AlignmentModel::new(
                                read_id,
                                &read_sequence,
                                &base_quality_scores,
                                &records
                            );
                            correct_rna_read(
                                *cluster_id,
                                read_name,
                                &alignment_model,
                                curr_cluster_graph_operations,
                                &chromosome_names_map,
                                gene_annotator,
                                options.max_correctable_event_len,
                                options.corrected_base_quality
                            )
                        })).unwrap_or_else(|_| {
                            log_info!("Failed to prepare read {}; emitting it unchanged.", read_name);
                            num_uncorrected.fetch_add(1, Ordering::Relaxed);
                            // A read stored without base qualities (SAM `*`) gets Q60 on every
                            // base, as AlignmentModel::new gives it.
                            let base_quality_scores: Vec<u8> = if base_quality_scores.is_empty() {
                                vec![60u8; read_sequence.len()]
                            } else {
                                base_quality_scores
                            };
                            CorrectedRNARead::new(
                                *cluster_id,
                                read_name.clone(),
                                read_sequence,
                                base_quality_scores
                            )
                        });

                        // Log progress.
                        let num_done: usize = num_corrected.fetch_add(1, Ordering::Relaxed) + 1;
                        if (num_done * 10) / num_reads != ((num_done - 1) * 10) / num_reads {
                            log_info!(
                                "Corrected {}% of reads ({}/{})",
                                (num_done * 100) / num_reads, num_done, num_reads
                            );
                        }

                        corrected_rna_read
                    }
                )
                .collect()
        });

        for corrected_read in corrected_reads.iter() {
            let quality: Vec<u8> = corrected_read.base_quality_scores.iter().map(|&q| q + b'!').collect();
            write_fastq_record(
                &mut writer,
                &corrected_read.read_name,
                corrected_read.sequence.as_bytes(),
                &quality
            );
        }
    }
    writer.finish().unwrap();
    std::fs::rename(&partial_fastq_file, output_fastq_file).unwrap();

    let num_uncorrected: usize = num_uncorrected.into_inner();
    if num_uncorrected > 0 {
        log_info!("{} of {} reads could not be corrected and were written as sequenced.", num_uncorrected, num_reads);
    }
    Ok(())
}


fn translate_cluster_operations(
    cluster_set: &RNAReadClusterSet,
    bam_chromosome_names_map: &BiMap<Box<str>, u16>
) -> HashMap<usize, HashSet<GraphOperation>> {
    let cluster_chromosome_names_map: &BiMap<Box<str>, u16> = cluster_set.get_chromosome_names_map();
    let translate = |chromosome_id: u16| -> Option<u16> {
        let name: &Box<str> = cluster_chromosome_names_map.get_by_right(&chromosome_id)?;
        bam_chromosome_names_map.get_by_left(name).copied()
    };
    cluster_set
        .get_clusters()
        .iter()
        .map(|cluster| {
            let operations: HashSet<GraphOperation> = cluster
                .get_variant_calls()
                .iter()
                .filter_map(|variant_call| {
                    let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
                    Some(GraphOperation::new(
                        translate(operation.get_chromosome_1())?,
                        operation.get_position_1(),
                        operation.get_strand_1().clone(),
                        operation.get_operation_type_1().clone(),
                        translate(operation.get_chromosome_2())?,
                        operation.get_position_2(),
                        operation.get_strand_2().clone(),
                        operation.get_operation_type_2().clone(),
                        operation.get_sequence().into(),
                        operation.get_variant_type().clone()
                    ))
                })
                .collect();
            (cluster.get_id(), operations)
        })
        .collect()
}


#[cfg(test)]
#[path = "../tests/pipeline/rna_read_correction.rs"]
mod tests;