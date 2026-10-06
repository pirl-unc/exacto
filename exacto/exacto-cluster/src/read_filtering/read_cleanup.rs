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
use exacto_core::prelude::FastaMap;
use std::collections::HashMap;

use crate::options::RNAVariantCallingOptions;
use crate::read_filtering::soft_clip::is_soft_clip_across_intron;
use crate::read_filtering::unplaced_tail::identify_unplaced_tail;


pub(crate) fn clean_transcript_model(
    transcript_model: &mut TranscriptModel,
    intron_boundary_index: &HashMap<(u16, u32), Vec<(u32, u32)>>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    calling: &RNAVariantCallingOptions
) {
    let read_length: u32 = transcript_model.get_alignment_model().num_bases();

    let unplaced_tail: Option<(u32, u32)> = if calling.unplaced_tail_removal {
        identify_unplaced_tail(
            transcript_model.get_variant_records(),
            transcript_model.get_splice_junctions(),
            read_length,
            calling.unplaced_tail_min_ins_len
        )
    } else {
        None
    };
    if let Some((tail_start, tail_end)) = unplaced_tail {
        transcript_model.retain_splice_junctions(|splice_junction| {
            splice_junction.read_position_2 < tail_start || splice_junction.read_position_1 > tail_end
        });
        transcript_model.retain_variant_records(|variant_record| {
            let read_position_1: u32 = variant_record.get_read_position_1().min(variant_record.get_read_position_2());
            let read_position_2: u32 = variant_record.get_read_position_1().max(variant_record.get_read_position_2());
            read_position_2 < tail_start || read_position_1 > tail_end
        });
    }

    let read_introns: Vec<(u16, u32, u32)> = transcript_model
        .get_splice_junctions_key()
        .iter()
        .map(|junction| junction.intron_span_key())
        .collect();
    transcript_model.retain_variant_records(|variant_record| {
        !calling.soft_clip_removal || !is_soft_clip_across_intron(
            variant_record,
            read_length,
            &read_introns,
            intron_boundary_index,
            chromosome_names_map,
            fasta_map,
            calling.soft_clip_max_boundary_distance,
            calling.soft_clip_bases_per_edit
        )
    });
}
