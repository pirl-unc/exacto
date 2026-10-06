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


use exacto_caller::prelude::*;


pub fn identify_unplaced_tail(
    variant_records: &[VariantRecord],
    splice_junctions: &[TranscriptModelSpliceJunction],
    read_length: u32,
    min_insertion_length: u32
) -> Option<(u32, u32)> {
    // Step 1. Get the read positions of every insertion.
    // Vec<(first read position, last read position)>
    let insertions: Vec<(u32, u32)> = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::Insertion)
        .map(|variant_record| (
            variant_record.get_read_position_1().min(variant_record.get_read_position_2()),
            variant_record.get_read_position_1().max(variant_record.get_read_position_2())
        ))
        .collect();

    // Step 2. Pair a terminal soft clip with a large insertion a short block away from it.
    for &(clip_start, clip_end) in insertions.iter() {
        let is_leading: bool = clip_start == 0;
        let is_trailing: bool = clip_end + 1 == read_length;
        if !is_leading && !is_trailing {
            continue;
        }
        let clip_length: u32 = clip_end - clip_start + 1;
        for &(insertion_start, insertion_end) in insertions.iter() {
            let insertion_length: u32 = insertion_end - insertion_start + 1;
            let is_terminal: bool = insertion_start == 0 || insertion_end + 1 == read_length;
            if is_terminal || insertion_length < min_insertion_length {
                continue;
            }
            let island_length: u32 = if is_trailing {
                clip_start.saturating_sub(insertion_end + 1)
            } else {
                insertion_start.saturating_sub(clip_end + 1)
            };
            if island_length == 0 || island_length >= insertion_length.min(clip_length) {
                continue;
            }

            // Step 3. Keep the pair when a splice junction lies between the insertion and the
            // block. A junction beside an insertion is stated over the aligned bases either side
            // of the inserted run.
            let is_spliced: bool = splice_junctions.iter().any(|splice_junction| {
                splice_junction.read_position_1 < insertion_start && insertion_end < splice_junction.read_position_2
            });
            if !is_spliced {
                continue;
            }
            return Some(if is_trailing {
                (insertion_start, clip_end)
            } else {
                (clip_start, insertion_end)
            });
        }
    }

    None
}


#[cfg(test)]
#[path = "../tests/read_filtering/unplaced_tail.rs"]
mod tests;