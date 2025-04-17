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


use std::any::Any;
use std::collections::{HashMap, HashSet};


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
/// use exacto_util::common::utilities::merge_regions;
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


pub fn count_non_overlapping_bases(
    a: &Vec<(Box<str>, u32, u32)>,
    b: &Vec<(Box<str>, u32, u32)>
) -> (u32,u32) {
    let mut a_bases: HashMap<Box<str>,HashSet<u32>> = HashMap::new();
    let mut b_bases: HashMap<Box<str>,HashSet<u32>> = HashMap::new();

    for (chromosome, start, end) in a.iter() {
        for i in *start..=*end {
            a_bases
                .entry(chromosome.clone())
                .or_insert_with(HashSet::new)
                .insert(i);
        }
    }

    for (chromosome, start, end) in b.iter() {
        for i in *start..=*end {
            b_bases
                .entry(chromosome.clone())
                .or_insert_with(HashSet::new)
                .insert(i);
        }
    }

    let mut a_only_bases: u32 = 0;
    for (chromosome, a_positions) in &a_bases {
        if b_bases.contains_key(chromosome) {
            let b_positions = b_bases.get(chromosome).unwrap();
            for i in a_positions {
                if b_positions.contains(i) == false {
                    a_only_bases += 1;
                }

            }
        } else {
            for _ in a_positions {
                a_only_bases += 1;
            }
        }
    }

    let mut b_only_bases: u32 = 0;
    for (chromosome, b_positions) in &b_bases {
        if a_bases.contains_key(chromosome) {
            let a_positions = a_bases.get(chromosome).unwrap();
            for i in b_positions {
                if a_positions.contains(i) == false {
                    b_only_bases += 1;
                }

            }
        } else {
            for _ in b_positions {
                b_only_bases += 1;
            }
        }
    }

    (a_only_bases, b_only_bases)
}