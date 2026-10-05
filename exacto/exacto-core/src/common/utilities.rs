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

pub fn calculate_l2_distance(vec1: &Vec<i8>, vec2: &Vec<i8>) -> f32 {
    assert_eq!(vec1.len(), vec2.len(), "Vectors must be the same length");
    let mut sum_sq_diff: f32 = 0.0;
    for (&a, &b) in vec1.iter().zip(vec2.iter()) {
        let diff: f32 = (a as f32) - (b as f32);
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

/// Count the bases covered by exactly one of the two interval sets.
///
/// Returns `(bases only in a, bases only in b)`. Endpoints are inclusive, a position is
/// counted once however many intervals cover it, and chromosomes are kept apart. An
/// interval whose start is past its end covers nothing.
///
/// Interval arithmetic rather than a coordinate bitmap. The bitmap version allocated one
/// bit per position from 0 to the highest end on the chromosome, so a locus a few megabases
/// into chr17 cost about a megabyte per side and twice that again for the clones the
/// AND/NOT needed — around 5 ms per call, which made this single function 86% of the RNA
/// read-modelling pass. Merging a few dozen exons gives the same answer from a sort.
pub fn count_non_overlapping_bases(
    a: &Vec<(Box<str>, u32, u32)>,
    b: &Vec<(Box<str>, u32, u32)>
) -> (u32, u32) {
    let a_merged: HashMap<&str, Vec<(u32, u32)>> = merge_regions_by_chromosome(a);
    let b_merged: HashMap<&str, Vec<(u32, u32)>> = merge_regions_by_chromosome(b);

    let count_private = |
        these: &HashMap<&str, Vec<(u32, u32)>>,
        those: &HashMap<&str, Vec<(u32, u32)>>
    | -> u32 {
        these
            .iter()
            .map(|(chromosome, intervals)| {
                let covered: u32 = intervals
                    .iter()
                    .map(|(start, end)| end - start + 1)
                    .sum();
                let shared: u32 = those
                    .get(chromosome)
                    .map_or(0, |other| count_intersecting_bases(intervals, other));
                covered - shared
            })
            .sum()
    };

    (count_private(&a_merged, &b_merged), count_private(&b_merged, &a_merged))
}

/// Group regions by chromosome and reduce each group to sorted, disjoint intervals.
///
/// Abutting intervals are merged along with overlapping ones. That does not change which
/// positions are covered, only how few pieces they are described in, which is what makes
/// the sweep in `count_intersecting_bases` linear.
fn merge_regions_by_chromosome(regions: &[(Box<str>, u32, u32)]) -> HashMap<&str, Vec<(u32, u32)>> {
    let mut by_chromosome: HashMap<&str, Vec<(u32, u32)>> = HashMap::new();
    for (chromosome, start, end) in regions.iter() {
        if start > end {
            continue;
        }
        by_chromosome
            .entry(chromosome.as_ref())
            .or_default()
            .push((*start, *end));
    }

    for intervals in by_chromosome.values_mut() {
        intervals.sort_unstable();

        let mut merged: Vec<(u32, u32)> = Vec::with_capacity(intervals.len());
        for &(start, end) in intervals.iter() {
            match merged.last_mut() {
                Some(last) if start <= last.1.saturating_add(1) => last.1 = last.1.max(end),
                _ => merged.push((start, end))
            }
        }
        *intervals = merged;
    }

    by_chromosome
}

/// Bases covered by both interval lists, each of which must be sorted and disjoint.
fn count_intersecting_bases(a: &[(u32, u32)], b: &[(u32, u32)]) -> u32 {
    let mut total: u32 = 0;
    let mut i: usize = 0;
    let mut j: usize = 0;

    while i < a.len() && j < b.len() {
        let start: u32 = a[i].0.max(b[j].0);
        let end: u32 = a[i].1.min(b[j].1);
        if start <= end {
            total += end - start + 1;
        }

        // Retire whichever interval ends first; the other may still meet the next one.
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }

    total
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

pub fn interval_contains(start_1: u32, end_1: u32, start_2: u32, end_2: u32) -> bool {
    start_1 <= start_2 && end_2 <= end_1
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
pub fn merge_regions(regions: Vec<(isize, isize)>) -> Vec<(isize, isize)> {
    let mut regions_: Vec<(isize,isize)> = regions.clone();
    regions_.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    if regions_.is_empty() {
        return Vec::new();
    }
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
