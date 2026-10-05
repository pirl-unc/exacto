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


use exacto_core::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::prelude::*;


pub struct NonsenseMediatedDecayPredictor<'a> {
    pub translation_strategy: &'a TranslationStrategy,
    pub start_codons: &'a HashSet<&'a str>,
    pub distance_threshold: u32
}

impl NonsenseMediatedDecayPredictor<'_> {
    pub fn predict(&self, transcript_model: &TranscriptModel) -> Vec<NonsenseMediatedDecayCall> {
        // Step 1. Get the last exon-exon junction read_position_1.
        // `read_position_1` is the last placed base of the exon 5' of it.
        // Read coordinates:  ... 198  199 | 200  201 ...
        //                          exon 1 | exon 2
        //                                 ^
        // read_position_1 = 199           |
        // read_position_2 = 200
        let last_junction: Option<ReadPosition> = transcript_model
            .get_splice_junctions()
            .iter()
            .map(|junction| junction.read_position_1)
            .max();

        // Step 2. Get the transcript sequence.
        let read_sequence: String = transcript_model.get_alignment_model().get_read_sequence();

        // Step 3. Identify open reading frames.
        let orfs: Vec<(Box<str>, ReadPosition, ReadPosition, u32)> = identify_open_reading_frames(
            &read_sequence,
            self.translation_strategy,
            self.start_codons
        );

        // Step 4. Classify nonsense mediated decay.
        orfs.iter()
            .map(|(_, orf_start, orf_end, _)| NonsenseMediatedDecayCall {
                orf_start: *orf_start,
                orf_end: *orf_end,
                verdict: classify_stop_codon(*orf_end, last_junction, self.distance_threshold)
            })
            .collect()
    }
}


fn classify_stop_codon(
    stop_codon_end: ReadPosition,
    last_junction: Option<ReadPosition>,
    nmd_distance_threshold: u32
) -> NonsenseMediatedDecayVerdict {
    use NonsenseMediatedDecayVerdict::*;

    // Step 1. Is there a junction downstream of the stop codon? 
    // Without one (an unspliced model, or a stop codon in the last exon) no exon junction complex survives termination.
    let Some(downstream_junction) = last_junction.filter(|&junction| junction >= stop_codon_end) else {
        return NotPredicted {
            distance_to_last_junction: None
        };
    };

    // Step 2. Is it far enough downstream that the terminating ribosome does not displace it?
    let distance: u32 = downstream_junction - stop_codon_end;
    if distance > nmd_distance_threshold {
        Predicted {
            distance_to_last_junction: distance
        }
    } else {
        NotPredicted {
            distance_to_last_junction: Some(distance)
        }
    }
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct NonsenseMediatedDecayCall {
    pub orf_start: ReadPosition,
    pub orf_end: ReadPosition,
    pub verdict: NonsenseMediatedDecayVerdict
}


#[cfg(test)]
#[path = "../tests/transcript/nonsense_mediated_decay_prediction.rs"]
mod tests;