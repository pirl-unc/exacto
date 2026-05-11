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
use exacto_core::prelude::*;
use rayon::prelude::*;
use statrs::distribution::{Binomial, DiscreteCDF};
use statrs::distribution::{Discrete, Hypergeometric};
use std::cmp::{min,max};
use std::collections::{BTreeMap,HashMap,HashSet};
use std::sync::Arc;

use crate::prelude::*;


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
/// * `stranded` is a boolean indicating whether the variants are stranded. If true, then the
/// variants must be on the same strand. If false, then the variants can be on the same or
/// opposite strands.
pub fn cluster_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    depths_map: &Arc<HashMap<Box<str>, Vec<u32>>>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    num_threads: usize,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    stranded: bool
) -> Vec<VariantCall> {
    // Step 1. Split variant records by chromosome
    // Key = (chromosome_1 ID, chromosome_2 ID)
    let mut variant_records_map: HashMap<(u16, u16, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> = split_variant_records(
        variant_records
            .into_iter()
            .map(|rc| Arc::new((*rc).clone()))
            .collect(),
        num_threads
    );

    // Step 2. Identify variant calls
    let mut variant_calls: Vec<VariantCall> = Vec::new();
    for curr_variant_records in variant_records_map.values() {
        // Identify local clusters of variant records using sweep line algorithm
        let clusters: Vec<VariantRecordCluster> = sweep_clusters(
            curr_variant_records,
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_intrachromosomal_distance_tau,
            max_intrachromosomal_distance,
            max_interchromosomal_distance,
            num_threads,
            stranded
        );

        // Create variant calls
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let curr_variant_calls: Vec<VariantCall> = thread_pool.install(|| {
            clusters
                .par_chunks(((clusters.len() + num_threads - 1) / num_threads).max(1))
                .flat_map(|curr_clusters| {
                    let mut curr_variant_calls: Vec<VariantCall> = Vec::new();
                    for curr_cluster in curr_clusters.iter() {
                        let mut curr_variant_call: VariantCall = VariantCall::new(0);
                        for variant_record in curr_cluster.variant_records.iter() {
                            curr_variant_call.add_variant_record((**variant_record).clone());
                        }
                        let variant_record: &VariantRecord = curr_variant_call.get_consensus_record().0;
                        let total_depth: u32 = get_variant_position_total_depth(
                            &depths_map,
                            chromosome_names_map.get_by_right(&variant_record.get_chromosome_1()).unwrap(),
                            variant_record.get_position_1(),
                            variant_record.get_operation_1(),
                            chromosome_names_map.get_by_right(&variant_record.get_chromosome_2()).unwrap(),
                            variant_record.get_position_2(),
                            variant_record.get_operation_2(),
                            variant_record.get_standardized_sequence().as_str()
                        );
                        curr_variant_call.set_total_depth(total_depth as i32);
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

pub fn compute_min_read_support(
    total_depth: u64,
    expected_variant_allele_fraction: f64,      // 0.25 for tumor, 0.5 for normal
    expected_mutation_rate: f64,                // 1e-6 for tumor, 1e-3 for normal
    expected_sequencing_error: f64,             // 0.01 for PacBio and Illumina, 0.05 for ONT (R9)
    max_f1_fraction: f64,                       // e.g. 0.99
    max_fpr: f64                                // e.g. 1e-6
) -> (u32, f64, f64, f64, f64) {
    assert!(total_depth > 0);
    assert!((0.0..=1.0).contains(&expected_variant_allele_fraction) && expected_variant_allele_fraction > 0.0);
    assert!((0.0..=1.0).contains(&expected_mutation_rate) && expected_mutation_rate > 0.0);
    assert!((0.0..=1.0).contains(&expected_sequencing_error));
    assert!((0.0..=1.0).contains(&max_f1_fraction));
    assert!((0.0..=1.0).contains(&max_fpr));

    let binom_f = Binomial::new(expected_variant_allele_fraction, total_depth).unwrap();
    let binom_e = Binomial::new(expected_sequencing_error, total_depth).unwrap();

    let mut f1_values: Vec<f64> = Vec::with_capacity(total_depth as usize);
    let mut recalls: Vec<f64> = Vec::with_capacity(total_depth as usize);
    let mut fpr_values: Vec<f64> = Vec::with_capacity(total_depth as usize);
    let mut precisions: Vec<f64> = Vec::with_capacity(total_depth as usize);

    for k in 1..=total_depth {
        let recall: f64 = 1.0 - binom_f.cdf(k - 1);
        let fpr: f64    = 1.0 - binom_e.cdf(k - 1);
        let precision: f64 = (recall * expected_mutation_rate) / ((recall * expected_mutation_rate) + (fpr * (1.0 - expected_mutation_rate)));
        let f1: f64 = f1_score(precision, recall);

        recalls.push(recall);
        fpr_values.push(fpr);
        precisions.push(precision);
        f1_values.push(f1);
    }

    // Eligible k: precision >= min_precision
    let eligible: Vec<usize> = (0..(total_depth as usize))
        .filter(|&idx| fpr_values[idx] <= max_fpr && f1_values[idx].is_finite())
        .collect();

    if eligible.is_empty() {
        return (u32::MAX, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    }

    // Peak F1 among eligible k
    let (peak_idx, f1_peak) = eligible
        .iter()
        .copied()
        .fold((eligible[0], f1_values[eligible[0]]), |best, i| {
            let v = f1_values[i];
            if v > best.1 { (i, v) } else { best }
        });

    let threshold: f64 = max_f1_fraction * f1_peak;

    // Among eligible ks, choose the smallest k whose F1 is within threshold of the eligible peak
    let best_idx = eligible
        .iter()
        .copied()
        .filter(|&i| f1_values[i] >= threshold)
        .min()
        .unwrap_or(peak_idx);

    let optimal_read_support = (best_idx + 1) as u32;
    (
        optimal_read_support,
        f1_values[best_idx],
        recalls[best_idx],
        fpr_values[best_idx],
        precisions[best_idx]
    )
}

pub fn compute_min_read_support_index(
    max_depth: u64,
    max_slippage_repeat_length: u32,
    expected_variant_allele_fraction: f64,  // 0.25 for tumor, 0.5 for normal
    expected_mutation_rate: f64,            // 1e-6 for tumor, 1e-3 for normal
    expected_sequencing_error: f64,         // 0.001 for PacBio and Illumina, 0.05 for ONT
    expected_slippage_rate: f64,            // same as expected_sequencing_error
    max_f1_fraction: f64,                   // typically 0.95
    max_fpr: f64                            // typically 1e-6
) -> Vec<Vec<u32>> {
    assert!(max_depth > 0);
    assert!(max_slippage_repeat_length > 1);
    assert!(expected_variant_allele_fraction > 0.0);
    assert!(expected_variant_allele_fraction < 1.0);
    assert!(expected_mutation_rate > 0.0);
    assert!(expected_mutation_rate < 1.0);
    assert!(expected_sequencing_error >= 0.0);
    assert!(expected_sequencing_error < 1.0);
    assert!(expected_slippage_rate >= 0.0);
    assert!(expected_slippage_rate < 1.0);
    assert!(max_f1_fraction > 0.0);
    assert!(max_f1_fraction <= 1.0);
    assert!(max_fpr > 0.0);
    assert!(max_fpr <= 1.0);

    let mut read_support_index: Vec<Vec<u32>> = Vec::new();

    // No repeat
    let mut read_suports: Vec<u32> = Vec::new();
    for c in 1..=max_depth {
        let (min_read_support, f1, recall, fpr, precision) = compute_min_read_support(
            c,
            expected_variant_allele_fraction,
            expected_mutation_rate,
            expected_sequencing_error,
            max_f1_fraction,
            max_fpr
        );
        read_suports.push(min_read_support);
    }
    read_support_index.push(read_suports);

    // No 1-bp homopolymer repeat
    read_support_index.push(vec![]);

    // Repeats (from 2-bp)
    for r in 2..=max_slippage_repeat_length {
        let mut read_suports: Vec<u32> = Vec::new();
        let p_slippage: f64 = 1.0 - (1.0 - expected_slippage_rate).powf(r as f64 - 1.0);
        for c in 1..=max_depth {
            let (min_read_support, f1, recall, fpr, precision) = compute_min_read_support(
                c,
                expected_variant_allele_fraction,
                expected_mutation_rate,
                p_slippage,
                max_f1_fraction,
                max_fpr
            );
            read_suports.push(min_read_support);
        }
        read_support_index.push(read_suports);
    }

    read_support_index
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
/// * `stranded` is a boolean indicating whether the variants are stranded. If true, then the
/// variants must be on the same strand. If false, then the variants can be on the same or
/// opposite strands.
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
    apply_infinite_sites_assumption: bool,
    stranded: bool
) -> Vec<Arc<VariantRecord>> {
    if b.is_empty() {
        return a;
    }

    // Step 1. Create indices
    let mut position_snv_map: HashSet<(u16, u32)> = HashSet::new();
    let mut position_1_map: HashMap<(u16, u32, VariantType),BTreeMap<u32, HashSet<Arc<GraphOperation>>>> = HashMap::new();
    let mut position_2_map: HashMap<(u16, u32, VariantType),BTreeMap<u32, HashSet<Arc<GraphOperation>>>> = HashMap::new();
    for variant_record in b.iter() {
        let variant_type: VariantType = variant_record.get_variant_type().clone();
        if variant_type == VariantType::SingleNucleotideVariant {
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
                .insert(Arc::new(variant_record.graph_operation.clone()));
            position_2_map
                .entry((chr2,zipcode2,variant_type.clone()))
                .or_insert_with(BTreeMap::new)
                .entry(pos2)
                .or_insert_with(HashSet::new)
                .insert(Arc::new(variant_record.graph_operation.clone()));
        }
    }

    // Step 2. Identify differences
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
                        let variant_type: VariantType = variant_a.get_variant_type().clone();
                        if variant_type == VariantType::SingleNucleotideVariant {
                            if position_snv_map.contains(&(variant_a.get_chromosome_1(),variant_a.get_position_1())) {
                                return None;
                            } else {
                                return Some(variant_a.clone());
                            }
                        } else {
                            // Search nearby variant records for position_1
                            let chr1: u16 = variant_a.get_chromosome_1();
                            let pos1: u32 = variant_a.get_position_1();
                            let zipcode1: u32 = pos1 / bin_size;
                            let min_zipcode1: u32 = if zipcode1 > 0 { zipcode1 - 1 } else { 0 };
                            for zipcode in min_zipcode1..=(zipcode1 + 1) {
                                if let Some(btree) = position_1_map.get(&(chr1,zipcode,variant_type.clone())) {
                                    let results = btree.range(pos1.saturating_sub(max_distance)..=pos1 + max_distance);
                                    for (_, sequence_operations) in results {
                                        for sequence_operation in sequence_operations.iter() {
                                            if !is_different(
                                                &variant_a.graph_operation,
                                                sequence_operation,
                                                min_size_proportion,
                                                max_ins_norm_edit_distance,
                                                max_intrachromosomal_distance_tau,
                                                max_intrachromosomal_distance,
                                                max_interchromosomal_distance,
                                                apply_infinite_sites_assumption,
                                                stranded
                                            ) {
                                                return None;
                                            }
                                        }
                                    }
                                }
                            }

                            // Search nearby variant records for position_2
                            let chr2: u16 = variant_a.get_chromosome_2();
                            let pos2: u32 = variant_a.get_position_2();
                            let zipcode2: u32 = pos2 / bin_size;
                            let min_zipcode2: u32 = if zipcode2 > 0 { zipcode2 - 1 } else { 0 };
                            for zipcode in min_zipcode2..=(zipcode2 + 1) {
                                if let Some(btree) = position_2_map.get(&(chr2,zipcode,variant_type.clone())) {
                                    let results = btree.range(pos2 - max_distance..=pos2 + max_distance);
                                    for (_, sequence_operations) in results {
                                        for sequence_operation in sequence_operations.iter() {
                                            if !is_different(
                                                &variant_a.graph_operation,
                                                sequence_operation,
                                                min_size_proportion,
                                                max_ins_norm_edit_distance,
                                                max_intrachromosomal_distance_tau,
                                                max_intrachromosomal_distance,
                                                max_interchromosomal_distance,
                                                apply_infinite_sites_assumption,
                                                stranded
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

pub fn get_variant_position_total_depth(
    depths_map: &HashMap<Box<str>, Vec<u32>>,
    chromosome_1: &str,
    position_1: u32,
    operation_1: &GraphOperationType,
    chromosome_2: &str,
    position_2: u32,
    operation_2: &GraphOperationType,
    sequence: &str
) -> u32 {
    // SNV
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) == 1 &&
        sequence.len() == 1 {
        let index: usize = (position_1 + 1 - 1) as usize;
        let depth: u32 = *depths_map.get(chromosome_1).unwrap().get(index).unwrap();
        return depth;
    }

    // MNV
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) == sequence.len() as u32 &&
        sequence.len() >= 2 {
        let start: usize = (position_1 + 1 - 1) as usize;
        let end: usize = (position_2 - 1 - 1) as usize;
        let depths: Vec<u32> = depths_map.get(chromosome_1).unwrap()[start..=end].to_vec();
        let max_depth: u32 = *depths.iter().max().unwrap();
        return max_depth;
    }

    // INS
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) == 0 &&
        sequence.len() >= 1 {
        let start: usize = (position_1 - 1) as usize;
        let end: usize = (position_2 - 1) as usize;
        let depths: Vec<u32> = depths_map.get(chromosome_1).unwrap()[start..=end].to_vec();
        let max_depth: u32 = *depths.iter().max().unwrap();
        return max_depth;
    }

    // DEL
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) >= 1  &&
        sequence.len() == 0 {
        let index_1: usize = (position_1 - 1) as usize;
        let index_2: usize = (position_2 - 1) as usize;
        let depth_1: u32 = *depths_map.get(chromosome_1).unwrap().get(index_1).unwrap();
        let depth_2: u32 = *depths_map.get(chromosome_1).unwrap().get(index_2).unwrap();
        if depth_1 < depth_2 {
            return depth_2;
        } else {
            return depth_1;
        }
    }

    // BND
    let index_1: usize = (position_1 - 1) as usize;
    let index_2: usize = (position_2 - 1) as usize;
    let depth_1: u32 = *depths_map.get(chromosome_1).unwrap().get(index_1).unwrap();
    let depth_2: u32 = *depths_map.get(chromosome_2).unwrap().get(index_2).unwrap();
    if depth_1 < depth_2 {
        depth_2
    } else {
        depth_1
    }
}

pub fn has_strand_bias(
    alt_fwd: u64,
    alt_rev: u64,
    ref_fwd: u64,
    ref_rev: u64,
    alpha: f64
) -> bool {
    if alt_fwd == 0 || alt_rev == 0 {
        return true;
    }

    let alt_total: u64 = alt_fwd + alt_rev;
    let ref_total: u64 = ref_fwd + ref_rev;

    // Fisher's exact test needs some REF evidence to compare against.
    if ref_total == 0 {
        return false;
    }

    let p: f64 = perform_fisher_exact_test_two_sided(
        alt_fwd,
        alt_rev,
        ref_fwd,
        ref_rev
    );
    
    p < alpha
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
/// * `stranded` is a boolean indicating whether the variants are stranded. If true, then the
/// variants must be on the same strand. If false, then the variants can be on the same or
/// opposite strands.
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
    max_interchromosomal_distance: u32,
    stranded: bool
) -> bool {
    // Records from the same read ID cannot be clustered
    let a_read_id: usize = a.get_read_id();
    let b_read_id: usize = b.get_read_id();
    if a_read_id == b_read_id {
        return false;
    }

    // Variant types must be the same
    let a_variant_type: &VariantType = a.get_variant_type();
    let b_variant_type: &VariantType = b.get_variant_type();
    if a_variant_type != b_variant_type {
        return false;
    }

    // Orientations must be the same
    let a_op1: &GraphOperationType = a.get_operation_1();
    let a_op2: &GraphOperationType = a.get_operation_2();
    let b_op1: &GraphOperationType = b.get_operation_1();
    let b_op2: &GraphOperationType = b.get_operation_2();
    if a_op1 != b_op1 || a_op2 != b_op2 {
        return false;
    }

    // Chromosomes 1 and 2 must be the same
    let a_chr1: u16 = a.get_chromosome_1();
    let a_chr2: u16 = a.get_chromosome_2();
    let b_chr1: u16 = b.get_chromosome_1();
    let b_chr2: u16 = b.get_chromosome_2();
    if a_chr1 != b_chr1 || a_chr2 != b_chr2 {
        return false;
    }

    // Strand pairing must be the same
    let a_s1: &Strand = a.get_strand_1();
    let a_s2: &Strand = a.get_strand_2();
    let b_s1: &Strand = b.get_strand_1();
    let b_s2: &Strand = b.get_strand_2();
    match (a_s1, a_s2) {
        (Strand::Forward, Strand::Forward) => {
            if b_s1 != b_s2 {
                return false;
            }
        },
        (Strand::Forward, Strand::Reverse) => {
            if b_s1 == b_s2 {
                return false;
            }
        },
        (Strand::Reverse, Strand::Forward) => {
            if b_s1 == b_s2 {
                return false;
            }
        },
        (Strand::Reverse, Strand::Reverse) => {
            if b_s1 != b_s2 {
                return false;
            }
        },
        _ => {
            panic!("Invalid strand combination: {} and {}", a.get_strand_1().as_str(), a.get_strand_2().as_str());
        }
    }

    // If stranded, strands must match exactly
    if stranded && (a_s1 != b_s1 || a_s2 != b_s2) {
        return false;
    }

    // The size proportion must be within the allowed limit
    let a_size: isize = a.get_variant_size();
    let b_size: isize = b.get_variant_size();
    let has_valid_sizes: bool = a_size >= 0 && b_size >= 0;
    if has_valid_sizes {
        let min_size: f32 = a_size.min(b_size) as f32;
        let max_size: f32 = a_size.max(b_size) as f32;
        if max_size > 0.0 {
            let size_proportion: f32 = min_size / max_size;
            if size_proportion < min_size_proportion {
                return false;
            }
        }
    }

    // The breakpoint distances must be close where proximity is a function of the variant size
    let max_distance: u32 = if !has_valid_sizes {
        max_interchromosomal_distance
    } else if a_chr1 == b_chr2 {
        let size_u32: u32 = a_size.max(b_size) as u32;
        if size_u32 == 1 {
            0
        } else {
            let tau = max_intrachromosomal_distance_tau as f32;
            let size_f = size_u32 as f32;
            let factor = 1.0_f32 - (-size_f / tau).exp();
            (max_intrachromosomal_distance as f32 * factor).ceil() as u32
        }
    } else {
        max_interchromosomal_distance
    };

    let a_pos1: u32 = a.get_position_1();
    let a_pos2: u32 = a.get_position_2();
    let b_pos1: u32 = b.get_position_1();
    let b_pos2: u32 = b.get_position_2();

    let distance_1: u32 = a_pos1.abs_diff(b_pos1);
    let distance_2: u32 = a_pos2.abs_diff(b_pos2);

    if distance_1 > max_distance || distance_2 > max_distance {
        return false;
    }

    // The normalized edit distance must be within the allowed limit if they are both insertions
    if *a_variant_type == VariantType::Insertion {
        let a_sequence: String = a.graph_operation.get_standardized_sequence();
        let b_sequence: String = b.graph_operation.get_standardized_sequence();
        let edit_distance: f32 = edit_distance(a_sequence.as_str(), b_sequence.as_str()) as f32;
        let max_len: f32 = a_sequence.len().max(b_sequence.len()) as f32;
        let normalized_edit_distance: f32 = edit_distance / max_len;
        if normalized_edit_distance > max_ins_norm_edit_distance {
            return false;
        }
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
/// * `stranded` is a boolean indicating whether the variants are stranded. If true, then the
/// variants must be on the same strand. If false, then the variants can be on the same or
/// opposite strands.
///
/// # Returns:
///
/// * True is `a` is different from `b`. False otherwise.
pub fn is_different(
    a: &GraphOperation,
    b: &GraphOperation,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    apply_infinite_sites_assumption: bool,
    stranded: bool
) -> bool {
    // Fast early exit for different chromosomes
    if (a.get_chromosome_1() != b.get_chromosome_1() && a.get_chromosome_2() != b.get_chromosome_2()) &&
        (a.get_chromosome_1() != b.get_chromosome_2() && a.get_chromosome_2() != b.get_chromosome_1()) {
        return true;
    }

    // Strand comparison
    if stranded {
        if a.get_strand_1() != b.get_strand_1() || a.get_strand_2() != b.get_strand_2() {
            return true;
        }
    }

    // Precompute reusable properties
    let size_a: isize = a.get_variant_size();
    let size_b: isize = b.get_variant_size();
    let max_size: isize = max(size_a, size_b);
    let min_size: isize = min(size_a, size_b);
    let size_proportion: f32 = if max_size > 0 {
        min_size as f32 / max_size as f32
    } else {
        1.0
    };

    // Infinite sites assumption
    if apply_infinite_sites_assumption {
        if (a.get_chromosome_1() == b.get_chromosome_1() && a.get_position_1() == b.get_position_1()) ||
            (a.get_chromosome_2() == b.get_chromosome_2() && a.get_position_2() == b.get_position_2()) {
            return false;
        }
    }

    let pos1_distance: u32 = a.get_position_1().abs_diff(b.get_position_1());
    let pos2_distance: u32 = a.get_position_2().abs_diff(b.get_position_2());

    // Translocation
    if a.get_variant_type().clone() == VariantType::Translocation {
        let cross_pos1_distance: u32 = a.get_position_1().abs_diff(b.get_position_2());
        let cross_pos2_distance: u32 = a.get_position_2().abs_diff(b.get_position_1());
        if (a.get_chromosome_1() == b.get_chromosome_1() && pos1_distance <= max_interchromosomal_distance) ||
            (a.get_chromosome_2() == b.get_chromosome_2() && pos2_distance <= max_interchromosomal_distance) ||
            (a.get_chromosome_1() == b.get_chromosome_2() && cross_pos1_distance <= max_interchromosomal_distance) ||
            (a.get_chromosome_2() == b.get_chromosome_1() && cross_pos2_distance <= max_interchromosomal_distance) {
            return false;
        } else {
            return true;
        }
    }

    let max_distance: u32 = calculate_max_distance(
        max(size_a, size_b) as f32, 
        max_intrachromosomal_distance_tau, 
        max_intrachromosomal_distance
    );

    // Insertion
    if a.get_variant_type().clone() == VariantType::Insertion {
        let a_sequence: String = a.get_standardized_sequence();
        let b_sequence: String = b.get_standardized_sequence();
        let edit_distance: f32 = edit_distance(a_sequence.as_str(), b_sequence.as_str()) as f32;
        let max_seq_length: f32 = max(a.get_sequence_length(), b.get_sequence_length()) as f32;
        let normalized_edit_distance: f32 = edit_distance / max_seq_length;
        if a.get_chromosome_1() == b.get_chromosome_1() &&
            pos1_distance <= max_distance &&
            size_proportion >= min_size_proportion &&
            normalized_edit_distance <= max_ins_norm_edit_distance {
            return false;
        } else {
            return true;
        }
    }

    // Deletion, breakpoint or MNV
    if (a.get_chromosome_1() == b.get_chromosome_1() && pos1_distance <= max_distance && size_proportion >= min_size_proportion) ||
        (a.get_chromosome_2() == b.get_chromosome_2() && pos2_distance <= max_distance && size_proportion >= min_size_proportion) {
        return false;
    }

    true
}

pub fn is_repeat_indel(
    variant_record: &VariantRecord,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> (bool, u32) {
    let variant_type: &VariantType = variant_record.get_variant_type();
    if variant_type != &VariantType::Insertion && variant_type != &VariantType::Deletion {
        return (false, 0);
    }

    let chromosome: &str = chromosome_names_map
        .get_by_right(&variant_record.get_chromosome_1())
        .unwrap();

    let position_1: u32 = variant_record.get_position_1();
    let position_2: u32 = variant_record.get_position_2();

    // Helper functions
    let get_fasta_sequence = |start: u32, end: u32| -> Box<str> {
        if start >= 1 && end <= fasta_map.get_length(chromosome) as u32 {
            fasta_map
                .get_sequence(chromosome, start as usize, end as usize)
                .to_uppercase()
                .into_boxed_str()
        } else {
            "".into()
        }
    };
    let count_left_homopolymer = |mut position: u32, base: &str| -> u32 {
        let mut n: u32 = 0;
        while position >= 1 {
            if &*get_fasta_sequence(position, position) == base {
                n += 1;
                if position == 1 {
                    break;
                }
                position -= 1;
            } else {
                break;
            }
        }
        n
    };
    let count_right_homopolymer = |mut position: u32, base: &str| -> u32 {
        let mut n: u32 = 0;
        loop {
            if &*get_fasta_sequence(position, position) == base {
                n += 1;
                position = position.saturating_add(1);
            } else {
                break;
            }
        }
        n
    };
    let count_left_dinucleotide_repeat = |mut position: u32, motif: &str| -> u32 {
        let mut n: u32 = 0;
        while position >= 2 {
            let s: Box<str> = get_fasta_sequence(position - 1, position);
            if s.len() != 2 {
                break;
            }
            if &*s == motif {
                n += 2;
                position -= 2;
            } else {
                break;
            }
        }
        n
    };
    let count_right_dinucleotide_repeat = |mut position: u32, motif: &str| -> u32 {
        let mut n: u32 = 0;
        loop {
            let s: Box<str> = get_fasta_sequence(position, position + 1);
            if s.len() != 2 {
                break;
            }
            if &*s == motif {
                n += 2;
                position = position.saturating_add(2);
            } else {
                break;
            }
        }
        n
    };

    // Insertion
    if variant_type == &VariantType::Insertion {
        let sequence: String = variant_record.graph_operation.get_standardized_sequence();
        assert!(sequence.len() > 0);

        // Homopolymer insertion
        if is_homopolymer_sequence(sequence.as_str()) {
            let base: &str = sequence.chars().next().map(|c| &sequence[..c.len_utf8()]).unwrap();

            // Count the length of the left homopolymer sequence, if any
            let left_homopolymer_length: u32 = count_left_homopolymer(position_1, base);

            // Count the length of the right homopolymer sequence, if any
            let right_homopolymer_length: u32 = count_right_homopolymer(position_2, base);

            if left_homopolymer_length > 0 || right_homopolymer_length > 0 {
                return (true, left_homopolymer_length + sequence.len() as u32 + right_homopolymer_length);
            }
        }

        // Dinucleotide repeat insertion
        if sequence.len() == 1 {
            // Count the length of the left dinucleotide repeat sequence, if any
            let left_sequence: Box<str> = get_fasta_sequence(position_1 - 1, position_1);
            let left_repeat_length: u32 = count_left_dinucleotide_repeat(position_1, &*left_sequence);

            // Count the length of the right dinucleotide repeat sequence, if any
            let right_sequence: Box<str> = get_fasta_sequence(position_2, position_2 + 1);
            let right_repeat_length: u32 = count_right_dinucleotide_repeat(position_2, &*right_sequence);

            if left_repeat_length > 0 || right_repeat_length > 0 {
                return (true, left_repeat_length + sequence.len() as u32 + right_repeat_length);
            }
        } else {
            let result: Option<String> = get_dinucleotide_repeat_motif(sequence.as_str());
            if result.is_some() {
                let motif: String = result.unwrap();
                let motif_reversed: String = motif.chars().rev().collect();

                // Count the length of the left dinucleotide repeat sequence, if any
                let left_repeat_length_1: u32 = count_left_dinucleotide_repeat(position_1, motif.as_str());
                let left_repeat_length_2: u32 = count_left_dinucleotide_repeat(position_1, motif_reversed.as_str());

                // Count the length of the right dinucleotide repeat sequence, if any
                let right_repeat_length_1: u32 = count_right_dinucleotide_repeat(position_2, motif.as_str());
                let right_repeat_length_2: u32 = count_right_dinucleotide_repeat(position_2, motif_reversed.as_str());

                let length_1: u32 = if left_repeat_length_1 > 0 || right_repeat_length_1 > 0 {
                    left_repeat_length_1 + sequence.len() as u32 + right_repeat_length_1
                } else {
                    0
                };

                let length_2: u32 = if left_repeat_length_2 > 0 || right_repeat_length_2 > 0 {
                    left_repeat_length_2 + sequence.len() as u32 + right_repeat_length_2
                } else {
                    0
                };

                if length_1 > 0 || length_2 > 0 {
                    return (true, length_1.max(length_2));
                }
            }
        }
    }

    // Deletion
    if variant_type == &VariantType::Deletion {
        let sequence: Box<str> = get_fasta_sequence(position_1 + 1, position_2 - 1);
        assert!(sequence.len() > 0);

        // Homopolymer deletion
        if is_homopolymer_sequence(&*sequence) {
            let base: &str = sequence.chars().next().map(|c| &sequence[..c.len_utf8()]).unwrap();

            // Count the length of the left homopolymer sequence, if any
            let left_homopolymer_length: u32 = count_left_homopolymer(position_1, base);

            // Count the length of the right homopolymer sequence, if any
            let right_homopolymer_length: u32 = count_right_homopolymer(position_2, base);

            if left_homopolymer_length > 0 || right_homopolymer_length > 0 {
                return (true, left_homopolymer_length + sequence.len() as u32 + right_homopolymer_length);
            }
        }

        // Dinucleotide repeat deletion
        if sequence.len() == 1 {
            // Count the length of the left dinucleotide repeat sequence, if any
            let left_sequence: Box<str> = get_fasta_sequence(position_1 - 1, position_1);
            let left_repeat_length: u32 = count_left_dinucleotide_repeat(position_1, &*left_sequence);

            // Count the length of the right dinucleotide repeat sequence, if any
            let right_sequence: Box<str> = get_fasta_sequence(position_2, position_2 + 1);
            let right_repeat_length: u32 = count_right_dinucleotide_repeat(position_2, &*right_sequence);

            if left_repeat_length > 0 || right_repeat_length > 0 {
                return (true, left_repeat_length + sequence.len() as u32 + right_repeat_length);
            }
        } else {
            let result: Option<String> = get_dinucleotide_repeat_motif(&*sequence);
            if result.is_some() {
                let motif: String = result.unwrap();
                let motif_reversed: String = motif.chars().rev().collect();

                // Count the length of the left dinucleotide repeat sequence, if any
                let left_repeat_length_1: u32 = count_left_dinucleotide_repeat(position_1, motif.as_str());
                let left_repeat_length_2: u32 = count_left_dinucleotide_repeat(position_1, motif_reversed.as_str());

                // Count the length of the right dinucleotide repeat sequence, if any
                let right_repeat_length_1: u32 = count_right_dinucleotide_repeat(position_2, motif.as_str());
                let right_repeat_length_2: u32 = count_right_dinucleotide_repeat(position_2, motif_reversed.as_str());

                let length_1: u32 = if left_repeat_length_1 > 0 || right_repeat_length_1 > 0 {
                    left_repeat_length_1 + sequence.len() as u32 + right_repeat_length_1
                } else {
                    0
                };

                let length_2: u32 = if left_repeat_length_2 > 0 || right_repeat_length_2 > 0 {
                    left_repeat_length_2 + sequence.len() as u32 + right_repeat_length_2
                } else {
                    0
                };

                if length_1 > 0 || length_2 > 0 {
                    return (true, length_1.max(length_2));
                }
            }
        }
    }

    (false, 0)
}


/// Two-sided Fisher exact test p-value for a 2x2 table.
///
/// Table:
///             FWD    REV
/// ALT          a      b
/// REF          c      d
///
/// Returns p-value in [0, 1].
pub fn perform_fisher_exact_test_two_sided(
    a: u64,
    b: u64,
    c: u64,
    d: u64
) -> f64 {
    let r1: u64 = a + b;    // ALT total
    let r2: u64 = c + d;    // REF total
    let c1: u64 = a + c;    // FWD total
    let c2: u64 = b + d;    // REV total
    let n: u64 = r1 + r2;   // total

    // Population: n
    // Population successes: c1 (FWD reads)
    // Draws: r1 (ALT reads)
    let hg: Hypergeometric = Hypergeometric::new(n, c1, r1).expect("Invalid hypergeometric parameters for Fisher test");
    let observed_p: f64 = hg.pmf(a);

    // Feasible range:
    // max(0, r1 - c2) <= x <= min(r1, c1)
    let min_x: u64 = r1.saturating_sub(n - c1);
    let max_x: u64 = r1.min(c1);

    // Two-sided definition: sum probabilities <= observed probability
    let mut p_two_sided: f64 = 0.0;
    for x in min_x..=max_x {
        let px: f64 = hg.pmf(x);
        if px <= observed_p + 1e-12 {
            p_two_sided += px;
        }
    }

    p_two_sided.min(1.0)
}

/// Split variant records by chromosome.
///
/// # Returns:
///
/// * A HashMap where the key is (chromosome_1,chromosome_2) and the value is a vector of
/// VariantRecord objects.
pub fn split_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    num_threads: usize
) -> HashMap<(u16, u16, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> {
    // Step 1. Split variant records by chromosome
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut variant_records_map: HashMap<(u16, u16, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> = thread_pool.install(|| {
        variant_records
            .par_iter()
            .map(|variant_record| {
                let key = (
                    variant_record.get_chromosome_1(),
                    variant_record.get_chromosome_2(),
                    variant_record.get_variant_type().clone(),
                    variant_record.get_operation_1().clone(),
                    variant_record.get_operation_2().clone()
                );
                (key, Arc::clone(variant_record))
            })
            .fold(
                || HashMap::new(),
                |mut acc, (key, variant_record)| {
                    acc.entry(key).or_insert_with(Vec::new).push(variant_record);
                    acc
                }
            )
            .reduce(
                || HashMap::new(),
                |mut map1, map2| {
                    for (key, mut vec) in map2 {
                        map1.entry(key).or_insert_with(Vec::new).append(&mut vec);
                    }
                    map1
                }
            )
    });

    // Step 2. Sort variant records by position_1
    thread_pool.install(|| {
        variant_records_map
            .iter_mut() // Borrow mutable references instead of consuming the map
            .for_each(|(_, records)| {
                records.sort_by(|variant_record_1, variant_record_2| {
                    variant_record_1.get_position_1().cmp(&variant_record_2.get_position_1())
                });
            });
    });

    variant_records_map
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
    variant_records: &Vec<Arc<VariantRecord>>,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    num_threads: usize,
    stranded: bool
) -> Vec<VariantRecordCluster> {
    let n: usize = variant_records.len();
    if n == 0 {
        return Vec::new();
    }

    // Step 1. Ensure all records are on the same chromosomes
    let chromosome_1: u16 = variant_records[0].get_chromosome_1();
    let chromosome_2: u16 = variant_records[0].get_chromosome_2();
    for vr in variant_records.iter() {
        assert!(
            vr.get_chromosome_1() == chromosome_1,
            "Supplied variant_records must have the same chromosome_1 value."
        );
        assert!(
            vr.get_chromosome_2() == chromosome_2,
            "Supplied variant_records must have the same chromosome_2 value."
        );
    }

    // Step 2. Identify clusterable pairs with a sliding window over sorted positions.
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let max_distance: u32 = max(max_intrachromosomal_distance, max_interchromosomal_distance);
    let positions: Vec<u32> = (0..n)
        .map(|i| variant_records[i].get_position_1())
        .collect();
    let mut j: usize = 0;
    let mut j_end: Vec<usize> = vec![n; n];
    for i in 0..n {
        if j < i + 1 {
            j = i + 1;
        }
        while j < n && positions[j].saturating_sub(positions[i]) <= max_distance {
            j += 1;
        }
        j_end[i] = j;
    }
    let indices: Vec<usize> = (0..variant_records.len()).collect();
    let pairs_list: Vec<(u32, Vec<u32>)> = thread_pool.install(|| {
        indices
            .par_iter()
            .map(|&i| {
                let mut mate_ids: Vec<u32> = Vec::new();
                let variant_record_i: &Arc<VariantRecord> = &variant_records[i];
                for j in (i + 1)..j_end[i] {
                    let variant_record_j: &Arc<VariantRecord> = &variant_records[j];
                    if is_clusterable(
                        variant_record_i.as_ref(),
                        variant_record_j.as_ref(),
                        min_size_proportion,
                        max_ins_norm_edit_distance,
                        max_intrachromosomal_distance_tau,
                        max_intrachromosomal_distance,
                        max_interchromosomal_distance,
                        stranded
                    ) {
                        mate_ids.push(j as u32);
                    }
                }
                (i as u32, mate_ids)
            })
            .collect()
    });

    // Step 4. Union-Find over indices
    let mut uf: UnionFind = UnionFind::new();
    for (i, mate_ids) in pairs_list.iter() {
        for j in mate_ids.iter() {
            uf.union(*i, *j);
        }
    }

    // Step 5. Build clusters from Union-Find sets
    let mut clusters: Vec<VariantRecordCluster> = Vec::new();
    let mut in_cluster: Vec<bool> = vec![false; n];
    for ids in uf.get_clusters().iter() {
        if ids.is_empty() {
            continue;
        }
        let mut cluster: VariantRecordCluster = VariantRecordCluster::new(
            chromosome_1,
            chromosome_2,
            u32::MAX,
            u32::MIN,
            u32::MAX,
            u32::MIN
        );

        for &idx in ids.iter() {
            let vr: Arc<VariantRecord> = Arc::clone(&variant_records[idx as usize]);

            // Update bounds
            let p1: u32 = vr.get_position_1();
            let p2: u32 = vr.get_position_2();
            if p1 < cluster.min_position_1 {
                cluster.min_position_1 = p1;
            }
            if p1 > cluster.max_position_1 {
                cluster.max_position_1 = p1;
            }
            if p2 < cluster.min_position_2 {
                cluster.min_position_2 = p2;
            }
            if p2 > cluster.max_position_2 {
                cluster.max_position_2 = p2;
            }

            cluster.variant_records.push(vr);
            in_cluster[idx as usize] = true;
        }

        clusters.push(cluster);
    }

    // Step 6. Add singleton clusters
    for (idx, vr) in variant_records.iter().enumerate() {
        if !in_cluster[idx] {
            let mut cluster = VariantRecordCluster::new(
                vr.get_chromosome_1(),
                vr.get_chromosome_2(),
                vr.get_position_1(),
                vr.get_position_1(),
                vr.get_position_2(),
                vr.get_position_2()
            );
            cluster.variant_records.push(Arc::clone(vr));
            clusters.push(cluster);
        }
    }

    clusters
}
