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


use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};
use serde::{Deserialize, Serialize};
use exacto_util::prelude::Transcript;


#[derive(Debug,Serialize,Deserialize)]
    pub struct ReferenceTranscriptMatch {
    pub reference_transcript: Transcript,
    pub num_overlap_bases: u32,
    pub num_transcript_only_bases: u32,
    pub num_reference_only_bases: u32,
    pub score: f32
}

impl PartialEq for ReferenceTranscriptMatch {
    fn eq(&self, other: &Self) -> bool {
        self.reference_transcript == other.reference_transcript &&
            self.num_overlap_bases == other.num_overlap_bases &&
            self.num_transcript_only_bases == other.num_transcript_only_bases &&
            self.num_reference_only_bases == other.num_reference_only_bases &&
            self.score == other.score
    }
}

impl Eq for ReferenceTranscriptMatch {}

impl Hash for ReferenceTranscriptMatch {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_transcript.hash(state);
        self.num_overlap_bases.hash(state);
        self.num_transcript_only_bases.hash(state);
        self.num_reference_only_bases.hash(state);
    }
}

impl ReferenceTranscriptMatch {
    pub fn new(
        reference_transcript: Transcript,
        num_overlap_bases: u32,
        num_transcript_only_bases: u32,
        num_reference_only_bases: u32,
        score: f32
    ) -> Self {
        Self {
            reference_transcript: reference_transcript,
            num_overlap_bases: num_overlap_bases,
            num_transcript_only_bases: num_transcript_only_bases,
            num_reference_only_bases: num_reference_only_bases,
            score: score
        }
    }
}

impl Clone for ReferenceTranscriptMatch {
    fn clone(&self) -> Self {
        ReferenceTranscriptMatch {
            reference_transcript: self.reference_transcript.clone(),
            num_overlap_bases: self.num_overlap_bases,
            num_transcript_only_bases: self.num_transcript_only_bases,
            num_reference_only_bases: self.num_reference_only_bases,
            score: self.score
        }
    }
}
