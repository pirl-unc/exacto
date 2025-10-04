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


use bitvec::prelude::*;
use std::any::Any;
use std::collections::HashMap;
use sysinfo::System;

use crate::log_info;


pub fn calculate_cosine_similarity(vec1: &Vec<i8>, vec2: &Vec<i8>) -> f64 {
    assert_eq!(vec1.len(), vec2.len(), "Vectors must be the same length");
    
    let (mut dot_product, mut norm_a_sq, mut norm_b_sq) = (0.0, 0.0, 0.0);
    for (&a, &b) in vec1.iter().zip(vec2.iter()) {
        let a_f64 = a as f64;
        let b_f64 = b as f64;
        dot_product += a_f64 * b_f64;
        norm_a_sq += a_f64 * a_f64;
        norm_b_sq += b_f64 * b_f64;
    }
    if norm_a_sq == 0.0 || norm_b_sq == 0.0 {
        0.0
    } else {
        dot_product / (norm_a_sq.sqrt() * norm_b_sq.sqrt())
    }
}

pub fn calculate_l2_distance(vec1: &Vec<i8>, vec2: &Vec<i8>) -> f64 {
    assert_eq!(vec1.len(), vec2.len(), "Vectors must be the same length");

    let mut sum_sq_diff = 0.0;
    for (&a, &b) in vec1.iter().zip(vec2.iter()) {
        let diff = (a - b) as f64;
        sum_sq_diff += diff * diff;
    }

    sum_sq_diff.sqrt()
}

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

/// Clone a Box<dyn Any>
pub fn clone_boxed_any(value: &Box<dyn Any>) -> Box<dyn Any> {
    if let Some(int_value) = value.downcast_ref::<i32>() {
        Box::new(int_value.clone()) as Box<dyn Any>
    } else if let Some(isize_value) = value.downcast_ref::<isize>() {
        Box::new(isize_value.clone()) as Box<dyn Any>
    } else if let Some(usize_value) = value.downcast_ref::<usize>() {
        Box::new(usize_value.clone()) as Box<dyn Any>
    } else if let Some(usize_value) = value.downcast_ref::<u8>() {
        Box::new(usize_value.clone()) as Box<dyn Any>
    } else if let Some(float_value) = value.downcast_ref::<f64>() {
        Box::new(float_value.clone()) as Box<dyn Any>
    } else if let Some(boxed_str_value) = value.downcast_ref::<Box<str>>() {
        Box::new(boxed_str_value.clone()) as Box<dyn Any>
    } else if let Some(bool_value) = value.downcast_ref::<bool>() {
        Box::new(bool_value.clone()) as Box<dyn Any>
    } else if let Some(vec) = value.downcast_ref::<Vec<Box<str>>>() {
        Box::new(vec.clone()) as Box<dyn Any>
    } else {
        panic!("Unsupported type for cloning in Box<dyn Any>");
    }
}

/// Find overlapping regions between two regions.
///
/// # Parameters:
///
/// * `segment_a` is a tuple of (start,end).
/// * `segment_b` is a tuple of (start,end).
///
/// # Returns:
///
/// * Option<(start,end)>.
/// * If there is no overlapping region between the two segments, returns `None`.
pub fn find_overlap(segment_a: (isize,isize), segment_b: (isize,isize)) -> Option<(isize,isize)> {
    let (a_start, a_end) = segment_a;
    let (b_start, b_end) = segment_b;
    let overlap_start = a_start.max(b_start);
    let overlap_end = a_end.min(b_end);
    if overlap_start <= overlap_end {
        Some((overlap_start, overlap_end))
    } else {
        None
    }
}

pub fn find_overlapping_region(intervals: &[(u32, u32)], value: u32) -> Option<(u32, u32)> {
    let is_sorted = intervals.is_sorted_by_key(|&(start, _)| start);
    assert!(is_sorted == true, "intervals must be sorted.");
    intervals.binary_search_by(|&(start, end)| {
        if value < start {
            std::cmp::Ordering::Greater
        } else if value > end {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    }).ok().map(|idx| intervals[idx])
}

pub fn is_overlapping(intervals: &[(u32, u32)], value: u32) -> bool {
    let is_sorted = intervals.is_sorted_by_key(|&(start, _)| start);
    assert!(is_sorted == true, "intervals must be sorted.");
    intervals.binary_search_by(|&(start, end)| {
        if value < start {
            std::cmp::Ordering::Greater
        } else if value > end {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    }).is_ok()
}

/// Check if two regions overlap.
///
/// Parameters:
///
/// * `start_1` is the start of region 1.
/// * `end_1` is the end of a region 1.
/// * `start_2` is the start of region 2.
/// * `end_2` is the end of region 2.
pub fn overlaps(start_1: isize, end_1: isize, start_2: isize, end_2: isize) -> bool {
    // De Morgan's law on checking for non-overlapping regions
    if start_1 <= end_2 && end_1 >= start_2 {
        true
    } else {
        false
    }
}

/// Merge a list of regions.
///
/// # Example
/// ```rust
/// use exacto_core::common::utilities::merge_regions;
///
/// let regions = vec![(1, 5), (2, 6), (8, 10), (9, 12)];
/// let merged = merge_regions(regions);
/// assert_eq!(merged, vec![(1, 6), (8, 12)]);
/// ```
pub fn merge_regions(regions: Vec<(isize, isize)>) -> Vec<(isize,isize)> {
    let mut regions_: Vec<(isize,isize)> = regions.clone();
    regions_.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut merged_regions = Vec::new();
    let mut current_region = regions_[0];
    for region in regions_.iter().skip(1) {
        if region.0 <= current_region.1 + 1 {
            // If the current region overlaps or is contiguous, extend it
            current_region.1 = current_region.1.max(region.1);
        } else {
            // Otherwise, finalize the current region and start a new one
            merged_regions.push(current_region);
            current_region = *region;
        }
    }
    merged_regions.push(current_region);
    merged_regions
}

pub fn count_common_bases(
    a: &Vec<(Box<str>, u32, u32)>,
    b: &Vec<(Box<str>, u32, u32)>
) -> u32 {
    let mut intervals: Vec<(u32, u32)> = Vec::new();

    // Step 1: Filter matching chromosomes and collect overlaps
    for (chr_a, start_a, end_a) in a {
        for (chr_b, start_b, end_b) in b {
            if chr_a == chr_b {
                let start = std::cmp::max(*start_a, *start_b);
                let end = std::cmp::min(*end_a, *end_b);
                if start <= end {
                    intervals.push((start, end));
                }
            }
        }
    }

    // Step 2: Merge overlapping intervals
    if intervals.is_empty() {
        return 0;
    }

    intervals.sort_by_key(|&(start, _)| start);
    let mut merged = vec![intervals[0]];

    for &(start, end) in &intervals[1..] {
        let last = merged.last_mut().unwrap();
        if start <= last.1 + 1 {
            last.1 = std::cmp::max(last.1, end);
        } else {
            merged.push((start, end));
        }
    }

    // Step 3: Compute inclusive base coverage
    merged.iter().map(|(start, end)| end - start + 1).sum()
}

pub fn count_union_bases(
    a: &Vec<(Box<str>, u32, u32)>,
    b: &Vec<(Box<str>, u32, u32)>
) -> u32 {
    let mut intervals: Vec<(Box<str>, u32, u32)> = Vec::new();
    intervals.extend_from_slice(a);
    intervals.extend_from_slice(b);

    // Group intervals by chromosome
    let mut unioned_len: u32 = 0;
    use std::collections::HashMap;
    let mut chromosome_intervals: HashMap<Box<str>,Vec<(u32, u32)>> = HashMap::new();

    for (chromosome, start, end) in intervals {
        chromosome_intervals
            .entry(chromosome.clone())
            .or_default()
            .push((start, end));
    }

    for (_chr, mut intervals) in chromosome_intervals {
        if intervals.is_empty() {
            continue;
        }

        // Sort and merge intervals for each chromosome
        intervals.sort_by_key(|&(start, _)| start);
        let mut merged: Vec<(u32, u32)> = vec![intervals[0]];

        for &(start, end) in &intervals[1..] {
            let last = merged.last_mut().unwrap();
            if start <= last.1 + 1 {
                last.1 = std::cmp::max(last.1, end);
            } else {
                merged.push((start, end));
            }
        }

        unioned_len += merged.iter().map(|(start, end)| end - start + 1).sum::<u32>();
    }

    unioned_len
}


// pub fn count_non_overlapping_bases(
//     a: &Vec<(Box<str>, u32, u32)>,
//     b: &Vec<(Box<str>, u32, u32)>
// ) -> (u32,u32) {
//     let mut a_bases: HashMap<&str, HashSet<u32>> = HashMap::new();
//     let mut b_bases: HashMap<&str, HashSet<u32>> = HashMap::new();
//
//     for (chromosome, start, end) in a.iter() {
//         a_bases
//             .entry(chromosome)
//             .or_insert_with(HashSet::new)
//             .extend(*start..=*end);
//     }
//     for (chromosome, start, end) in b.iter() {
//         b_bases
//             .entry(chromosome)
//             .or_insert_with(HashSet::new)
//             .extend(*start..=*end);
//     }
//
//     let mut a_only_bases: u32 = 0;
//     for (chromosome, a_set) in &a_bases {
//         if let Some(b_set) = b_bases.get(chromosome) {
//             a_only_bases += a_set.difference(b_set).count() as u32;
//         } else {
//             a_only_bases += a_set.len() as u32;
//         }
//     }
//
//     let mut b_only_bases = 0;
//     for (chromosome, b_set) in &b_bases {
//         if let Some(a_set) = a_bases.get(chromosome) {
//             b_only_bases += b_set.difference(a_set).count() as u32;
//         } else {
//             b_only_bases += b_set.len() as u32;
//         }
//     }
//
//     (a_only_bases, b_only_bases)
// }

// pub fn count_non_overlapping_bases(
//     a: &Vec<(Box<str>, u32, u32)>,
//     b: &Vec<(Box<str>, u32, u32)>
// ) -> (u32,u32) {
//     let mut a_bases: HashMap<&str, HashSet<u32>> = HashMap::new();
//     let mut b_bases: HashMap<&str, HashSet<u32>> = HashMap::new();
// 
//     for (chromosome, start, end) in a.iter() {
//         a_bases
//             .entry(chromosome)
//             .or_insert_with(HashSet::new)
//             .extend(*start..=*end);
//     }
//     for (chromosome, start, end) in b.iter() {
//         b_bases
//             .entry(chromosome)
//             .or_insert_with(HashSet::new)
//             .extend(*start..=*end);
//     }
// 
//     let mut a_only_bases: u32 = 0;
//     for (chromosome, a_set) in &a_bases {
//         if let Some(b_set) = b_bases.get(chromosome) {
//             a_only_bases += a_set.difference(b_set).count() as u32;
//         } else {
//             a_only_bases += a_set.len() as u32;
//         }
//     }
// 
//     let mut b_only_bases = 0;
//     for (chromosome, b_set) in &b_bases {
//         if let Some(a_set) = a_bases.get(chromosome) {
//             b_only_bases += b_set.difference(a_set).count() as u32;
//         } else {
//             b_only_bases += b_set.len() as u32;
//         }
//     }
// 
//     (a_only_bases, b_only_bases)
// }

pub fn count_non_overlapping_bases(
    a: &Vec<(Box<str>, u32, u32)>,
    b: &Vec<(Box<str>, u32, u32)>
) -> (u32, u32) {
    // Step 1. Determine max end per chromosome
    let mut chromosome_max_end: HashMap<&str, u32> = HashMap::new();
    for (chromosome, start, end) in a.iter().chain(b.iter()) {
        let entry = chromosome_max_end
            .entry(chromosome.as_ref())
            .or_insert(0);
        *entry = (*entry).max(*end);
    }

    // Step 2. Allocate bitvectors
    let mut a_map: HashMap<&str, BitVec> = HashMap::new();
    let mut b_map: HashMap<&str, BitVec> = HashMap::new();
    for (chromosome, start, end) in a.iter() {
        let length = (chromosome_max_end[chromosome.as_ref()] + 1) as usize;
        let bv = a_map.entry(chromosome).or_insert_with(|| bitvec![0; length]);
        for i in *start as usize..=*end as usize {
            bv.set(i, true);
        }
    }
    for (chromosome, start, end) in b.iter() {
        let length = (chromosome_max_end[chromosome.as_ref()] + 1) as usize;
        let bv = b_map.entry(chromosome).or_insert_with(|| bitvec![0; length]);
        for i in *start as usize..=*end as usize {
            bv.set(i, true);
        }
    }

    // Step 3. Count unique bases
    let mut a_only = 0u32;
    for (chromosome, a_bits) in &a_map {
        match b_map.get(chromosome) {
            Some(b_bits) => {
                a_only += (a_bits.clone() & !b_bits.clone()).count_ones() as u32;
            },
            None => {
                a_only += a_bits.count_ones() as u32;
            }
        }
    }

    let mut b_only = 0u32;
    for (chromosome, b_bits) in &b_map {
        match a_map.get(chromosome) {
            Some(a_bits) => {
                b_only += (b_bits.clone() & !a_bits.clone()).count_ones() as u32;
            },
            None => {
                b_only += b_bits.count_ones() as u32;
            }
        }
    }

    (a_only, b_only)
}
