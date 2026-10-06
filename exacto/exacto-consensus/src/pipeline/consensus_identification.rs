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


use abpoa_rs::AlignmentMode;
use exacto_core::prelude::*;
use exacto_core::log_info;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::prelude::*;


pub fn sort_clusters_by_num_reads(
    clusters: &HashMap<usize, HashSet<Box<str>>>
) -> Vec<(usize, &HashSet<Box<str>>)> {
    let mut sorted_clusters: Vec<(usize, &HashSet<Box<str>>)> = clusters
        .iter()
        .map(|(&cluster_id, read_names)| (cluster_id, read_names))
        .collect();
    sorted_clusters.sort_unstable_by_key(|&(cluster_id, read_names)| (Reverse(read_names.len()), cluster_id));
    sorted_clusters
}


pub fn identify_consensus_sequences(
    clusters: &HashMap<usize, HashSet<Box<str>>>,
    fastq_file: &str,
    options: &IdentifyConsensusSequencesOptions,
    num_threads: usize
) -> ConsensusSequenceSet {
    assert!(options.poa_gap_open_score >= 0, "poa_gap_open_score must be at least 0.");
    assert!(options.poa_gap_extend_score >= 1, "poa_gap_extend_score must be at least 1.");
    assert!(options.orientation_kmer_size <= 32, "orientation_kmer_size must be at most 32.");

    // Step 1. Build an index of the FASTQ file, and check that it holds every read of the clusters.
    log_info!("Building an offset index of the FASTQ file.");
    let offsets: HashMap<Box<str>, FastqOffset> = build_fastq_offset_index(fastq_file);
    let mut missing_read_names: Vec<&Box<str>> = clusters
        .values()
        .flatten()
        .filter(|read_name| !offsets.contains_key(*read_name))
        .collect();
    assert_eq!(missing_read_names.is_empty(), true);

    // Step 2. Sort the clusters by their number of reads, largest first.
    let sorted_clusters: Vec<(usize, &HashSet<Box<str>>)> = sort_clusters_by_num_reads(clusters);

    // Step 3. Identify the consensus sequence in each cluster.
    log_info!("Identify consensus sequence for each RNA cluster.");
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let total_clusters: usize = clusters.len();
    let total_reads: usize = clusters.values().map(|read_names| read_names.len()).sum();
    let processed_clusters: AtomicUsize = AtomicUsize::new(0);
    let processed_reads: AtomicUsize = AtomicUsize::new(0);
    let reversed_reads: AtomicUsize = AtomicUsize::new(0);
    let reoriented_clusters: AtomicUsize = AtomicUsize::new(0);
    let subsampled_clusters: AtomicUsize = AtomicUsize::new(0);
    // HashMap<cluster ID, consensus sequence>
    let consensus_by_cluster: HashMap<usize, Box<str>> = thread_pool.install(|| {
        sorted_clusters
            .into_iter()
            .par_bridge()
            .map_init(
                || FastqReader::open(fastq_file),   // one reader per worker thread, reused
                |reader, (cluster_id, read_names)| {
                    // Fetch the read sequences.
                    let aligned_read_names: Vec<&str> = subsample_read_names(
                        read_names,
                        options.max_reads_per_cluster,
                        (options.seed as u64).wrapping_add(cluster_id as u64)
                    );
                    if aligned_read_names.len() < read_names.len() {
                        subsampled_clusters.fetch_add(1, Ordering::Relaxed);
                    }
                    let mut sequences: Vec<Box<str>> = Vec::new();
                    for read_name in aligned_read_names {
                        let offset: FastqOffset = offsets[read_name];
                        let sequence: Box<str> = reader.get_read_sequence(read_name, offset);
                        sequences.push(sequence);
                    }
                    let num_reversed_reads: usize = orient_reads(&mut sequences, options.orientation_kmer_size);
                    if num_reversed_reads > 0 {
                        reversed_reads.fetch_add(num_reversed_reads, Ordering::Relaxed);
                        reoriented_clusters.fetch_add(1, Ordering::Relaxed);
                    }
                    sequences.sort_unstable_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
                    
                    // Perform partial order alignment.
                    let consensus: Box<str> = perform_partial_order_alignment(
                        &sequences,
                        AlignmentMode::Global,
                        options.poa_match_score,
                        options.poa_mismatch_score,
                        options.poa_gap_open_score,
                        options.poa_gap_extend_score
                    );

                    // Log the progress.
                    let clusters_done: usize = processed_clusters.fetch_add(1, Ordering::Relaxed) + 1;
                    let reads_before: usize = processed_reads.fetch_add(read_names.len(), Ordering::Relaxed);
                    let reads_done: usize = reads_before + read_names.len();
                    if reads_done * 10 / total_reads != reads_before * 10 / total_reads {
                        log_info!(
                            "Consensus progress: {}% of reads ({}/{}); {}/{} clusters done.",
                            reads_done * 100 / total_reads,
                            reads_done,
                            total_reads,
                            clusters_done,
                            total_clusters
                        );
                    }

                    (cluster_id, consensus)
                }
            )
            .collect()
    });
    if options.orientation_kmer_size > 0 {
        log_info!(
            "Reverse-complemented {} reads in {} clusters to the orientation of most reads of their cluster.",
            reversed_reads.load(Ordering::Relaxed),
            reoriented_clusters.load(Ordering::Relaxed)
        );
    }
    if options.max_reads_per_cluster > 0 {
        log_info!(
            "{} clusters of more than {} reads were aligned on {} of their reads, drawn at random.",
            subsampled_clusters.load(Ordering::Relaxed),
            options.max_reads_per_cluster,
            options.max_reads_per_cluster
        );
    }

    // Step 4. Merge the set, attaching each cluster's read names.
    log_info!("Merging the consensus sequences into a set.");
    let mut set: ConsensusSequenceSet = ConsensusSequenceSet::new();
    let mut cluster_ids: Vec<&usize> = clusters.keys().collect();
    cluster_ids.sort();
    for &cluster_id in cluster_ids {
        let read_names: &HashSet<Box<str>> = &clusters[&cluster_id];
        let consensus_sequence: Box<str> = consensus_by_cluster
            .get(&cluster_id)
            .cloned()
            .unwrap_or_else(|| "".into());
        set.add(ConsensusSequence::new(
            cluster_id,
            consensus_sequence,
            read_names.clone()
        ));
    }

    set
}


#[cfg(test)]
#[path = "../tests/pipeline/consensus_identification.rs"]
mod tests;