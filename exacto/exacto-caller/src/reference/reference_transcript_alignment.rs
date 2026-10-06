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


use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::{Alignment, AlignmentOperation};
use exacto_core::prelude::{reverse_complement, Strand};
use exacto_core::prelude::ReferenceTranscriptPosition;
use std::ops::Range;

use crate::prelude::ReferenceTranscriptSequence;


#[derive(Debug)]
pub(crate) struct ReferenceTranscriptPlacement {
    pub query_range: Range<u32>,                                // in the query as aligned (already reverse-complemented if !is_forward)
    pub reference_range: Range<ReferenceTranscriptPosition>,    // [start, end) reference transcript base positions
    pub is_forward: bool,                                       // query in reference transcript orientation, or its reverse complement
    pub score: i32
}


pub(crate) fn place_on_reference_transcript(
    query: &str,
    rts: &ReferenceTranscriptSequence,
    window: Range<ReferenceTranscriptPosition>,
    gap_open: i32,
    gap_extend: i32,
    k: u32,
    band_width: u32,
    min_score_fraction: f64,
    min_query_coverage: f64
) -> Option<ReferenceTranscriptPlacement> {
    if query.is_empty() || window.is_empty() {
        return None;
    }

    let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let sequence: String = rts.get_sequence().to_ascii_uppercase();
    let reference: &[u8] = &sequence.as_bytes()[window.start as usize..window.end as usize];
    let query_rc: Box<str> = reverse_complement(query);

    let mut aligner: Aligner<_> = Aligner::new(
        gap_open,
        gap_extend,
        score_fn,
        k as usize,
        band_width as usize
    );

    let mut best: Option<ReferenceTranscriptPlacement> = None;

    // Get the better alignment between the forward and reverse complement query.
    for (q, is_forward) in [(query, true), (&*query_rc, false)] {
        let alignment: Alignment = aligner.local(q.as_bytes(), reference);
        let aligned: usize = alignment.xend.saturating_sub(alignment.xstart);
        if aligned == 0
            || (aligned as f64) < min_query_coverage * query.len() as f64
            || (alignment.score as f64) < min_score_fraction * aligned as f64 {
            continue;
        }
        let candidate = ReferenceTranscriptPlacement {
            query_range: alignment.xstart as u32..alignment.xend as u32,
            reference_range: (window.start + alignment.ystart as u32)..(window.start + alignment.yend as u32),
            is_forward,
            score: alignment.score
        };
        if best.as_ref().map_or(true, |b| candidate.score > b.score) {
            best = Some(candidate);
        }
    }

    best
}


#[cfg(test)]
#[path = "../tests/reference/reference_transcript_alignment.rs"]
mod tests;