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

use crate::io::records::StitchedTranscriptRecord;
use crate::stitching::end_resolution::EndResolution;
use crate::stitching::stitched_transcript::StitchedTranscript;
use crate::stitching::stitched_transcript_set::StitchedTranscriptSet;


pub fn build_stitched_transcript_records<'a>(
    stitched_transcript_set: &'a StitchedTranscriptSet,
    read_names_map: &'a BiMap<Box<str>, usize>
) -> impl Iterator<Item = StitchedTranscriptRecord> + 'a {
    stitched_transcript_set
        .transcripts
        .iter()
        .map(move |transcript| build_stitched_transcript_record(transcript, read_names_map))
}


fn build_stitched_transcript_record(
    transcript: &StitchedTranscript,
    read_names_map: &BiMap<Box<str>, usize>
) -> StitchedTranscriptRecord {
    let read_name: &Box<str> = read_names_map
        .get_by_right(&transcript.read_id)
        .unwrap_or_else(|| panic!("read ID {} has no read name", transcript.read_id));

    let get_reference_ids = |resolution: &Option<EndResolution>| -> (Box<str>, Box<str>) {
        match resolution {
            Some(resolution) => (
                resolution.reference_transcript_match.get_reference_gene_id().into(),
                resolution.reference_transcript_match.get_reference_transcript_id().into()
            ),
            None => ("".into(), "".into())
        }
    };
    let (five_prime_reference_gene_id, five_prime_reference_transcript_id) =
        get_reference_ids(&transcript.five_prime_end_resolution);
    let (three_prime_reference_gene_id, three_prime_reference_transcript_id) =
        get_reference_ids(&transcript.three_prime_end_resolution);

    let five_prime_stitch_end: u32 = transcript
        .five_prime_stitch_end_index
        .unwrap_or(0) as u32;
    let three_prime_stitch_start: u32 = transcript
        .three_prime_stitch_start_index
        .unwrap_or(transcript.stitched_length as usize) as u32;

    let retained_read_start: u32 = match &transcript.five_prime_end_resolution {
        Some(resolution) => resolution.read_join_position,
        None => 0
    };
    let retained_read_end: u32 = match &transcript.three_prime_end_resolution {
        Some(resolution) => resolution.read_join_position + 1,
        None => transcript.original_length
    };

    StitchedTranscriptRecord {
        read_name: read_name.clone(),
        stitched_sequence: transcript.stitched_sequence.clone(),
        read_length: transcript.original_length,
        stitched_length: transcript.stitched_length,
        retained_read_start: retained_read_start,
        retained_read_end: retained_read_end,

        is_five_prime_stitched: transcript.is_five_prime_stitched(),
        five_prime_reference_gene_id: five_prime_reference_gene_id,
        five_prime_reference_transcript_id: five_prime_reference_transcript_id,
        five_prime_stitch_start: 0,
        five_prime_stitch_end: five_prime_stitch_end,
        five_prime_stitch_length: five_prime_stitch_end,

        is_three_prime_stitched: transcript.is_three_prime_stitched(),
        three_prime_reference_gene_id: three_prime_reference_gene_id,
        three_prime_reference_transcript_id: three_prime_reference_transcript_id,
        three_prime_stitch_start: three_prime_stitch_start,
        three_prime_stitch_end: transcript.stitched_length,
        three_prime_stitch_length: transcript.stitched_length - three_prime_stitch_start,

        unmodeled_reason: transcript.unmodeled_reason.unwrap_or("").into()
    }
}
