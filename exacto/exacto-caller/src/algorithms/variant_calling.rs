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
use edit_distance::edit_distance;
use exacto_util::prelude::*;
use interavl::IntervalTree;
use rayon::prelude::*;
use std::cmp::{min,max};
use std::collections::{BTreeMap,HashMap,HashSet};
use std::sync::Arc;
use sysinfo::System;

use crate::common::constants::*;
use crate::log_info;
use crate::structs::sequence_operation::SequenceOperation;
use crate::structs::variant_call::VariantCall;
use crate::structs::variant_record::VariantRecord;
use crate::structs::variant_record_cluster::VariantRecordCluster;


pub fn capture_memory_usage(message: &str) {
    let mut sys = System::new_all();
    sys.refresh_all();
    let pid = sysinfo::get_current_pid().unwrap();
    if let Some(process) = sys.process(pid) {
        let memory_usage = process.memory();
        let memory_usage_gb = memory_usage as f64 / (1024.0 * 1024.0 * 1024.0);
        log_info!("{}: {:.2} GB", message, memory_usage_gb);
    } else {
        log_info!("Could not get process memory usage");
    }
}

fn calculate_max_distance(size: f32, tau: u32, max_distance: u32) -> u32 {
    (max_distance as f32 * (1.0 - f32::exp(-size / tau as f32))).ceil() as u32
}

/// Cluster variant records.
///
/// # Parameters:
///
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_tau` is the `tau` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_intrachromosomal_distance` is the `d_max` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
pub fn cluster_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    num_threads: usize,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32
) -> Vec<VariantCall> {
    // Step 1. Split variant records
    // Key = (chromosome_1 ID, chromosome ID)
    let mut variant_records_map: HashMap<(u16,u16),Vec<Arc<VariantRecord>>> = split_variant_records_by_chromosome(
        variant_records
            .into_iter()
            .map(|rc| Arc::new((*rc).clone()))
            .collect(),
        num_threads
    );

    // Step 2. Sort variant records by position_1
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        variant_records_map
            .iter_mut() // Borrow mutable references instead of consuming the map
            .for_each(|(_, records)| {
                records.sort_by(|variant_record_1, variant_record_2| {
                    variant_record_1.get_position_1().cmp(&variant_record_2.get_position_1())
                });
            });
    });

    // Step 3. Identify variant calls
    let mut variant_calls: Vec<VariantCall> = Vec::new();
    for ((chromosome_1,chromosome_2), curr_variant_records) in variant_records_map.iter() {
        // Identify local clusters of variant records using sweep line algorithm
        let clusters: Vec<VariantRecordCluster> = sweep_clusters(
            curr_variant_records.clone(),
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_intrachromosomal_distance_tau,
            max_intrachromosomal_distance,
            max_interchromosomal_distance,
            num_threads
        );

        // Create variant calls
        let curr_variant_calls: Vec<VariantCall> = thread_pool.install(|| {
            clusters
                .par_chunks((clusters.len() + num_threads - 1) / num_threads)
                .flat_map(|curr_clusters| {
                    let mut curr_variant_calls: Vec<VariantCall> = Vec::new();
                    for curr_cluster in curr_clusters.iter() {
                        let mut curr_variant_call: VariantCall = VariantCall::new();
                        for variant_record in curr_cluster.variant_records.iter() {
                            curr_variant_call.add_variant_record((**variant_record).clone());
                        }
                        curr_variant_calls.push(curr_variant_call);
                    }
                    curr_variant_calls.into_par_iter()
                })
                .collect()
        });

        // Append variant calls
        variant_calls.extend(curr_variant_calls);
    }

    variant_calls
}

/// Diff b from a.
///
/// # Parameters:
///
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_tau` is the `tau` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_intrachromosomal_distance` is the `d_max` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
/// * If `apply_infinite_sites_assumption` is true, any `a` variant record that shares a breakpoint
/// (either position_1 or position_2) with any of the `b` variant record will be filtered out.
///
/// # Returns:
///
/// * A vector of VariantRecord objects specific to `a`.
pub fn diff_variant_records(
    a: Vec<Arc<VariantRecord>>,
    b: Vec<Arc<VariantRecord>>,
    bin_size: u32,
    num_threads: usize,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    apply_infinite_sites_assumption: bool
) -> Vec<Arc<VariantRecord>> {
    if b.is_empty() {
        return a;
    }

    // Step 1. Create indices
    log_info!("\tCreating indices");
    capture_memory_usage("\tBefore creating indices");
    let mut position_snv_map: HashSet<(u16,u32)> = HashSet::new();
    let mut position_1_map: HashMap<(u16,u32,SequenceOperationVariantTypes),BTreeMap<u32,HashSet<Arc<SequenceOperation>>>> = HashMap::new();
    let mut position_2_map: HashMap<(u16,u32,SequenceOperationVariantTypes),BTreeMap<u32,HashSet<Arc<SequenceOperation>>>> = HashMap::new();
    for variant_record in b.iter() {
        let variant_type: SequenceOperationVariantTypes = variant_record.get_variant_type();
        if variant_type == SequenceOperationVariantTypes::SingleNucleotideVariant {
            position_snv_map.insert((variant_record.get_chromosome_1(), variant_record.get_position_1()));
        } else {
            let chr1: u16 = variant_record.get_chromosome_1();
            let chr2: u16 = variant_record.get_chromosome_2();
            let pos1: u32 = variant_record.get_position_1();
            let pos2: u32 = variant_record.get_position_2();
            let zipcode1: u32 = pos1 / bin_size;
            let zipcode2: u32 = pos2 / bin_size;
            position_1_map
                .entry((chr1,zipcode1,variant_type.clone()))
                .or_insert_with(BTreeMap::new)
                .entry(pos1)
                .or_insert_with(HashSet::new)
                .insert(Arc::new(variant_record.sequence_operation.clone()));
            position_2_map
                .entry((chr2,zipcode2,variant_type.clone()))
                .or_insert_with(BTreeMap::new)
                .entry(pos2)
                .or_insert_with(HashSet::new)
                .insert(Arc::new(variant_record.sequence_operation.clone()));
        }
    }
    capture_memory_usage("\tAfter creating indices");
    log_info!("\t{} SNVs", position_snv_map.len());

    // Step 2. Identify differences
    log_info!("\tDiffing against interval trees");
    let max_distance: u32 = max(max_intrachromosomal_distance, max_interchromosomal_distance);
    let chunk_size = (a.len() + num_threads - 1) / num_threads;
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        a.par_chunks(chunk_size)
            .flat_map(|chunk| {
                chunk
                    .par_iter()
                    .filter_map(|variant_a| {
                        let variant_type: SequenceOperationVariantTypes = variant_a.get_variant_type();
                        if variant_type == SequenceOperationVariantTypes::SingleNucleotideVariant {
                            if position_snv_map.contains(&(variant_a.get_chromosome_1(),variant_a.get_position_1())) {
                                return None;
                            } else {
                                return Some(variant_a.clone());
                            }
                        } else {
                            // Search nearby variant records for position_1
                            let chr1 = variant_a.get_chromosome_1();
                            let pos1 = variant_a.get_position_1();
                            let zipcode1 = pos1 / bin_size;
                            let min_zipcode1 = if zipcode1 > 0 { zipcode1 - 1 } else { 0 };
                            for zipcode in min_zipcode1..=(zipcode1 + 1) {
                                if let Some(btree) = position_1_map.get(&(chr1,zipcode,variant_type.clone())) {
                                    let results = btree.range(pos1 - max_distance..=pos1 + max_distance);
                                    for (_, graph_operations) in results {
                                        for graph_operation in graph_operations.iter() {
                                            if !is_different(
                                                &variant_a.sequence_operation,
                                                graph_operation,
                                                min_size_proportion,
                                                max_ins_norm_edit_distance,
                                                max_intrachromosomal_distance_tau,
                                                max_intrachromosomal_distance,
                                                max_interchromosomal_distance,
                                                apply_infinite_sites_assumption
                                            ) {
                                                return None;
                                            }
                                        }
                                    }
                                }
                            }

                            // Search nearby variant records for position_2
                            let chr2 = variant_a.get_chromosome_2();
                            let pos2 = variant_a.get_position_2();
                            let zipcode2 = pos2 / bin_size;
                            let min_zipcode2 = if zipcode2 > 0 { zipcode2 - 1 } else { 0 };
                            for zipcode in min_zipcode2..=(zipcode2 + 1) {
                                if let Some(btree) = position_2_map.get(&(chr2,zipcode,variant_type.clone())) {
                                    let results = btree.range(pos2 - max_distance..=pos2 + max_distance);
                                    for (_, graph_operations) in results {
                                        for graph_operation in graph_operations.iter() {
                                            if !is_different(
                                                &variant_a.sequence_operation,
                                                graph_operation,
                                                min_size_proportion,
                                                max_ins_norm_edit_distance,
                                                max_intrachromosomal_distance_tau,
                                                max_intrachromosomal_distance,
                                                max_interchromosomal_distance,
                                                apply_infinite_sites_assumption
                                            ) {
                                                return None;
                                            }
                                        }
                                    }
                                }
                            }

                            Some(variant_a.clone())
                        }
                    })
            })
            .collect()
    })
}

/// Checks whether two VariantRecord objects can be clustered together.
///
/// # Parameters:
///
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_tau` is the `tau` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_intrachromosomal_distance` is the `d_max` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
///
/// # Returns:
///
/// True if the two variant records can be clustered. False otherwise.
pub fn is_clusterable(
    a: &VariantRecord,
    b: &VariantRecord,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32
) -> bool {
    // Chromosomes 1 and 2 must be the same
    if a.get_chromosome_1() != b.get_chromosome_1() || a.get_chromosome_2() != b.get_chromosome_2() {
        return false;
    }

    // Variant types must be the same
    if a.get_variant_type() != b.get_variant_type() {
        return false;
    }

    // Orientations must be the same
    if a.get_operation_1() != b.get_operation_1() || a.get_operation_2() != b.get_operation_2() {
        return false;
    }

    // The normalized edit distance must be within the allowed limit if they are both insertions
    if a.get_variant_type() == SequenceOperationVariantTypes::Insertion {
        let edit_distance = edit_distance(a.get_sequence(), b.get_sequence()) as f32;
        let max_size = f32::max(a.get_sequence_length() as f32, b.get_sequence_length() as f32);
        let normalized_edit_distance: f32 = edit_distance / max_size;
        if normalized_edit_distance > max_ins_norm_edit_distance {
            return false;
        }
    }

    // The breakpoint distances must be close where proximity is a function of the variant size
    let variant_size = f32::max(a.get_variant_size() as f32, b.get_variant_size() as f32);
    let max_distance: u32;
    if a.get_chromosome_1() == b.get_chromosome_2() {
        if variant_size == 1f32 {
            max_distance = 0;
        } else {
            max_distance = (max_intrachromosomal_distance as f32 * (1f32 - f32::exp(-1f32 * (variant_size / max_intrachromosomal_distance_tau as f32)))).ceil() as u32;
        }
    } else {
        max_distance = max_interchromosomal_distance;
    }
    let distance_1: u32 = (a.get_position_1() as f32 - b.get_position_1() as f32).abs() as u32;
    let distance_2: u32 = (a.get_position_2() as f32 - b.get_position_2() as f32).abs() as u32;
    if distance_1 > max_distance || distance_2 > max_distance {
        return false;
    }

    // The size proportion must be within the allowed limit
    let min_size = f32::min(a.get_variant_size() as f32, b.get_variant_size() as f32);
    let max_size = f32::max(a.get_variant_size() as f32, b.get_variant_size() as f32);
    let size_proportion: f32 = min_size / max_size;
    if size_proportion < min_size_proportion {
        return false;
    }

    // If none of the above conditions was met, then the two variant records can be clustered
    true
}

/// Check if a VariantRecord is different from b.
///
/// # Parameters:
///
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_tau` is the `tau` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_intrachromosomal_distance` is the `d_max` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
/// * If `apply_infinite_sites_assumption` is true, then if `a` shares a breakpoint
/// (either position_1 or position_2) with `b`, then the return value will be false.
///
/// # Returns:
///
/// * True is `a` is different from `b`. False otherwise.
pub fn is_different(
    graph_operation_1: &SequenceOperation,
    graph_operation_2: &SequenceOperation,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    apply_infinite_sites_assumption: bool
) -> bool {
    // Fast early exit for different chromosomes
    if graph_operation_1.chromosome_1 != graph_operation_2.chromosome_1 &&
        graph_operation_1.chromosome_2 != graph_operation_2.chromosome_2 {
        return true;
    }

    // Precompute reusable properties
    let size_a = graph_operation_1.get_variant_size();
    let size_b = graph_operation_2.get_variant_size();
    let max_size = max(size_a, size_b);
    let min_size = min(size_a, size_b);
    let size_proportion = min_size as f32 / max_size as f32;

    // Infinite sites assumption
    if apply_infinite_sites_assumption {
        if (graph_operation_1.chromosome_1 == graph_operation_2.chromosome_1 &&
            graph_operation_1.position_1 == graph_operation_2.position_1) ||
            (graph_operation_1.chromosome_2 == graph_operation_2.chromosome_2 &&
                graph_operation_1.position_2 == graph_operation_2.position_2) {
            return false;
        }
    }

    let pos1_distance = graph_operation_1.position_1.abs_diff(graph_operation_2.position_1);
    let pos2_distance = graph_operation_1.position_2.abs_diff(graph_operation_2.position_2);

    // Translocation
    if graph_operation_1.variant_type == SequenceOperationVariantTypes::Translocation {
        if (graph_operation_1.chromosome_1 == graph_operation_2.chromosome_1 && pos1_distance <= max_interchromosomal_distance) ||
            (graph_operation_1.chromosome_2 == graph_operation_2.chromosome_2 && pos2_distance <= max_interchromosomal_distance) {
            return false;
        } else {
            return true;
        }
    }

    let max_distance: u32 = calculate_max_distance(max(size_a, size_b) as f32, max_intrachromosomal_distance_tau, max_intrachromosomal_distance);

    // Insertion
    if graph_operation_1.variant_type == SequenceOperationVariantTypes::Insertion {
        let edit_distance: f32 = edit_distance(&*graph_operation_1.sequence, &*graph_operation_2.sequence) as f32;
        let max_seq_length: f32 = max(graph_operation_1.get_sequence_length(), graph_operation_2.get_sequence_length()) as f32;
        let normalized_edit_distance = edit_distance / max_seq_length;
        if graph_operation_1.chromosome_1 == graph_operation_2.chromosome_1 &&
            pos1_distance <= max_distance &&
            size_proportion >= min_size_proportion &&
            normalized_edit_distance <= max_ins_norm_edit_distance {
            return false;
        } else {
            return true;
        }
    }

    // Deletion, breakpoint or MNV
    if (graph_operation_1.chromosome_1 == graph_operation_2.chromosome_1 &&
        pos1_distance <= max_distance &&
        size_proportion >= min_size_proportion) ||
        (graph_operation_1.chromosome_2 == graph_operation_2.chromosome_2 &&
            pos2_distance <= max_distance &&
            size_proportion >= min_size_proportion) {
        return false;
    }

    true
}

/// Split variant records by chromosome.
///
/// # Returns:
///
/// * A HashMap where the key is (chromosome_1,chromosome_2) and the value is a vector of
/// VariantRecord objects.
pub fn split_variant_records_by_chromosome(
    variant_records: Vec<Arc<VariantRecord>>,
    num_threads: usize,
) -> HashMap<(u16, u16),Vec<Arc<VariantRecord>>> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        variant_records
            .par_iter()
            .map(|variant_record| {
                let key = (
                    variant_record.get_chromosome_1(),
                    variant_record.get_chromosome_2(),
                );
                (key, Arc::clone(variant_record))
            })
            .fold(
                || HashMap::new(),
                |mut acc, (key, variant_record)| {
                    acc.entry(key).or_insert_with(Vec::new).push(variant_record);
                    acc
                },
            )
            .reduce(
                || HashMap::new(),
                |mut map1, map2| {
                    for (key, mut vec) in map2 {
                        map1.entry(key).or_insert_with(Vec::new).append(&mut vec);
                    }
                    map1
                },
            )
    })
}

/// Sweep clusters.
///
/// # Parameters:
///
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_tau` is the `tau` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_intrachromosomal_distance` is the `d_max` term of the formula for determining the
/// intrachromosomal `max_breakpoint_distance (d)`:
/// `d = d_max * (1 - e^(-1*variant_size / tau))`. `variant_size` is the bigger of the two variant sizes.
/// Both breakpoint pairs must be within `max_breakpoint_distance (or d)` to be considered clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
pub fn sweep_clusters(
    variant_records: Vec<Arc<VariantRecord>>,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    num_threads: usize
) -> Vec<VariantRecordCluster> {
    // Step 1. Build a binary search tree
    let mut bst: BTreeMap<usize,Arc<VariantRecord>> = BTreeMap::new();
    for variant_record in variant_records.iter() {
        bst.insert(variant_record.get_position_1() as usize, Arc::clone(variant_record));
    }

    // Step 2. Ensure all records are on the same chromosomes
    let chromosome_1: u16 = variant_records[0].get_chromosome_1();
    let chromosome_2: u16 = variant_records[0].get_chromosome_2();
    for variant_record in variant_records.iter() {
        assert!(variant_record.get_chromosome_1() == chromosome_1, "Supplied variant_records must have the same chromosome_1 value.");
        assert!(variant_record.get_chromosome_2() == chromosome_2, "Supplied variant_records must have the same chromosome_2 value.");
    }

    // Step 3. Assign a unique ID to each VariantRecord
    let mut variant_records_map: BiMap<Arc<VariantRecord>,usize> = BiMap::new();
    for (id, variant_record) in variant_records.iter().enumerate() {
        variant_records_map.insert(Arc::clone(variant_record), id);
    }

    // Step 4. Identify clusterable pairs
    let max_distance: u32 = max(max_intrachromosomal_distance, max_interchromosomal_distance);
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let pairs: HashSet<(usize,usize)> = thread_pool.install(|| {
        variant_records
            .par_iter()
            .flat_map(|variant_record_1| {
                let min_position: usize = variant_record_1.get_position_1().saturating_sub(max_distance) as usize;
                let max_position: usize = (variant_record_1.get_position_1() + max_distance) as usize;
                let variant_records_: Vec<Arc<VariantRecord>> = bst.range(min_position..=max_position)
                    .map(|(_, variant_record)| Arc::clone(variant_record))
                    .collect();
                let mut curr_pairs: HashSet<(usize,usize)> = HashSet::new();
                for variant_record_2 in variant_records_.iter() {
                    if is_clusterable(
                        variant_record_1,
                        variant_record_2,
                        min_size_proportion,
                        max_ins_norm_edit_distance,
                        max_intrachromosomal_distance_tau,
                        max_intrachromosomal_distance,
                        max_interchromosomal_distance,
                    ) {
                        let variant_record_1_id: usize = *variant_records_map.get_by_left(variant_record_1).unwrap();
                        let variant_record_2_id: usize = *variant_records_map.get_by_left(variant_record_2).unwrap();
                        if variant_record_1_id < variant_record_2_id {
                            curr_pairs.insert((variant_record_1_id, variant_record_2_id));
                        } else {
                            curr_pairs.insert((variant_record_2_id, variant_record_1_id));
                        }
                    }
                }
                curr_pairs
            })
            .collect()
    });

    // Step 5. Identify clusters
    let mut uf: UnionFind = UnionFind::new();
    for pair in pairs.iter() {
        uf.union(pair.0, pair.1);
    }

    // Step 6. Get clusters from UnionFind
    let mut clusters: Vec<VariantRecordCluster> = Vec::new();
    let mut clustered_variant_records: HashSet<Arc<VariantRecord>> = HashSet::new();
    for variant_record_ids in uf.get_clusters() {
        let mut curr_variant_records: HashSet<Arc<VariantRecord>> = HashSet::new();
        let mut min_position_1: usize = usize::MAX;
        let mut min_position_2: usize = usize::MAX;
        let mut max_position_1: usize = usize::MIN;
        let mut max_position_2: usize = usize::MIN;
        for variant_record_id in variant_record_ids.iter() {
            let variant_record: &Arc<VariantRecord> = variant_records_map.get_by_right(variant_record_id).unwrap();
            curr_variant_records.insert(Arc::clone(variant_record));
            if (variant_record.get_position_1() as usize) < min_position_1 {
                min_position_1 = variant_record.get_position_1() as usize;
            }
            if (variant_record.get_position_1() as usize) > max_position_1 {
                max_position_1 = variant_record.get_position_1() as usize;
            }
            if (variant_record.get_position_2() as usize) < min_position_2 {
                min_position_2 = variant_record.get_position_2() as usize;
            }
            if (variant_record.get_position_2() as usize) > max_position_2 {
                max_position_2 = variant_record.get_position_2() as usize;
            }
        }
        let mut cluster = VariantRecordCluster::new(
            chromosome_1,
            chromosome_2,
            min_position_1 as u32,
            max_position_1 as u32,
            min_position_2 as u32,
            max_position_2 as u32,
        );
        for variant_record in curr_variant_records.iter() {
            cluster.add_variant_record(Arc::clone(variant_record));
            clustered_variant_records.insert(Arc::clone(variant_record));
        }
        clusters.push(cluster);
    }

    // Step 7. If a variant record was not included in any cluster, include them
    // as an independent VariantRecordCluster
    for variant_record in variant_records.iter() {
        if clustered_variant_records.contains(variant_record) == false {
            let mut cluster: VariantRecordCluster = VariantRecordCluster::new(
                variant_record.get_chromosome_1(),
                variant_record.get_chromosome_2(),
                variant_record.get_position_1(),
                variant_record.get_position_1(),
                variant_record.get_position_2(),
                variant_record.get_position_2()
            );
            cluster.add_variant_record(Arc::clone(variant_record));
            clusters.push(cluster);
        }
    }

    clusters
}
