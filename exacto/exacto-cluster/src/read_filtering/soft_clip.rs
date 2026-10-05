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
use exacto_core::prelude::{reverse_complement, FastaMap, Strand};
use std::collections::{BTreeSet, HashMap};


pub fn build_intron_boundary_index(
    summaries: &[RNAReadCharacterizationSummary]
) -> HashMap<(u16, u32), Vec<(u32, u32)>> {
    let introns: BTreeSet<(u16, u32, u32)> = summaries
        .iter()
        .flat_map(|summary| summary.splice_junctions.iter())
        .filter(|junction| junction.chromosome_1 == junction.chromosome_2)
        .map(|junction| junction.intron_span_key())
        .collect();
    let mut intron_boundary_index: HashMap<(u16, u32), Vec<(u32, u32)>> = HashMap::new();
    for (chromosome, intron_start, intron_end) in introns {
        intron_boundary_index.entry((chromosome, intron_start)).or_default().push((intron_start, intron_end));
        intron_boundary_index.entry((chromosome, intron_end)).or_default().push((intron_start, intron_end));
    }
    intron_boundary_index
}


pub fn is_soft_clip_across_intron(
    variant_record: &VariantRecord,
    read_length: u32,
    read_introns: &[(u16, u32, u32)],
    intron_boundary_index: &HashMap<(u16, u32), Vec<(u32, u32)>>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    max_boundary_distance: u32,
    num_bases_per_edit: u32
) -> bool {
    // Step 1. Keep the insertions made of a terminal soft clip.
    if *variant_record.get_variant_type() != VariantType::Insertion {
        return false;
    }
    let is_leading: bool = variant_record.get_read_position_1().min(variant_record.get_read_position_2()) == 0;
    let is_trailing: bool = variant_record.get_read_position_1().max(variant_record.get_read_position_2()) + 1 == read_length;
    if !is_leading && !is_trailing {
        return false;
    }
    let chromosome_id: u16 = variant_record.get_chromosome_1();
    let Some(chromosome) = chromosome_names_map.get_by_right(&chromosome_id) else {
        return false;
    };

    // Step 2. Spell the clip along the reference and find the side of the alignment it hangs off.
    let is_reverse: bool = *variant_record.get_strand_1() == Strand::Reverse;
    let clip: String = if is_reverse {
        reverse_complement(variant_record.get_sequence()).to_uppercase()
    } else {
        variant_record.get_sequence().to_uppercase()
    };
    let is_above: bool = is_leading == is_reverse;
    let anchor: i64 = if is_above {
        variant_record.get_position_1().min(variant_record.get_position_2()) as i64
    } else {
        variant_record.get_position_1().max(variant_record.get_position_2()) as i64
    };
    let max_edits: usize = if num_bases_per_edit == 0 {
        0
    } else {
        clip.len() / num_bases_per_edit as usize
    };
    let reach: i64 = (clip.len() + max_edits) as i64;
    let reference = |start: i64, end: i64| -> String {
        if start > end || end < 1 {
            return String::new();
        }
        fasta_map
            .try_get_sequence(chromosome, start.max(1) as usize, end as usize)
            .unwrap_or("")
            .to_uppercase()
    };

    // Step 3. Compare the clip to the reference across every intron with a boundary nearby.
    for boundary in (anchor - max_boundary_distance as i64)..=(anchor + max_boundary_distance as i64) {
        // The exonic base next to the intron is `boundary`; the intron starts one above it when
        // the clip hangs off the upper end of the alignment and ends one below it otherwise.
        let intron_base: i64 = if is_above {
            boundary + 1
        } else {
            boundary - 1
        };
        if intron_base < 1 {
            continue;
        }
        let Some(introns) = intron_boundary_index.get(&(chromosome_id, intron_base as u32)) else {
            continue;
        };
        for &(intron_start, intron_end) in introns.iter() {
            let (intron_start, intron_end): (i64, i64) = (intron_start as i64, intron_end as i64);
            if (is_above && intron_start != intron_base) || (!is_above && intron_end != intron_base) {
                continue;
            }

            // The exonic bases the read may go on from: the one across this intron, and the
            // one across every intron the read splices.
            let mut far_bases: BTreeSet<i64> = read_introns
                .iter()
                .filter(|(read_chromosome_id, _, _)| *read_chromosome_id == chromosome_id)
                .map(|&(_, read_intron_start, read_intron_end)| {
                    if is_above {
                        read_intron_end as i64 + 1
                    } else {
                        read_intron_start as i64 - 1
                    }
                })
                .collect();
            far_bases.insert(if is_above { intron_end + 1 } else { intron_start - 1 });

            for far_base in far_bases.into_iter() {
                let target: String = if is_above {
                    if anchor > boundary {
                        // The alignment ran into the intron: those bases have to be the first
                        // bases of the exon the read goes on from too.
                        let num_overrun: i64 = anchor - boundary;
                        if reference(intron_start, anchor) != reference(far_base, far_base + num_overrun - 1) {
                            continue;
                        }
                        reference(far_base + num_overrun, far_base + num_overrun + reach - 1)
                    } else {
                        reference(anchor + 1, boundary) + &reference(far_base, far_base + reach - 1)
                    }
                } else {
                    // Read away from the alignment, so the clip and the reference both start at
                    // the end of the alignment.
                    if anchor < boundary {
                        let num_overrun: i64 = boundary - anchor;
                        if reference(anchor, intron_end) != reference(far_base - num_overrun + 1, far_base) {
                            continue;
                        }
                        reference(far_base - num_overrun - reach + 1, far_base - num_overrun)
                            .chars()
                            .rev()
                            .collect()
                    } else {
                        (reference(far_base - reach + 1, far_base) + &reference(boundary, anchor - 1))
                            .chars()
                            .rev()
                            .collect()
                    }
                };
                let clip: String = if is_above {
                    clip.clone()
                } else {
                    clip.chars().rev().collect()
                };
                let min_length: usize = clip.len() - max_edits;
                let max_length: usize = (clip.len() + max_edits).min(target.len());
                if is_within_edits_of_a_prefix(clip.as_bytes(), target.as_bytes(), min_length, max_length, max_edits) {
                    return true;
                }
            }
        }
    }

    false
}


pub(crate) fn is_within_edits_of_a_prefix(clip: &[u8], target: &[u8], min_length: usize, max_length: usize, max_edits: usize) -> bool {
    if min_length > max_length {
        return false;
    }
    // Any distance beyond `max_edits` is held as `max_edits + 1`.
    let beyond: usize = max_edits + 1;
    // previous[j]: distance of the clip's first i - 1 bases to the target's first j.
    let mut previous: Vec<usize> = (0..=max_length).map(|j| j.min(beyond)).collect();
    let mut current: Vec<usize> = vec![beyond; max_length + 1];
    let (mut low, mut high): (usize, usize) = (0, max_edits.min(max_length));
    for i in 1..=clip.len() {
        (low, high) = (i.saturating_sub(max_edits), (i + max_edits).min(max_length));
        if low > high {
            return false;
        }
        if low == 0 {
            current[0] = i.min(beyond);
        } else {
            current[low - 1] = beyond;
        }
        for j in low.max(1)..=high {
            let substitution: usize = previous[j - 1] + usize::from(clip[i - 1] != target[j - 1]);
            current[j] = substitution.min(previous[j] + 1).min(current[j - 1] + 1).min(beyond);
        }
        if high < max_length {
            current[high + 1] = beyond;
        }
        std::mem::swap(&mut previous, &mut current);
    }
    (min_length.max(low)..=max_length.min(high)).any(|length| previous[length] <= max_edits)
}


#[cfg(test)]
#[path = "../tests/read_filtering/soft_clip.rs"]
mod tests;