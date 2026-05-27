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
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptModelSet {
    pub assembled_transcripts: HashSet<AssembledTranscript>,
    pub transcript_models: Vec<TranscriptModel>,
    pub read_names_map: BiMap<Box<str>, usize>,
    pub chromosome_names_map: BiMap<Box<str>, u16>
}

impl TranscriptModelSet {
    pub fn new() -> Self {
        Self {
            assembled_transcripts: HashSet::new(),
            transcript_models: Vec::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add_transcript_model(
        &mut self,
        assembled_transcript: AssembledTranscript,
        transcript_models: Vec<TranscriptModel>
    ) {
        self.assembled_transcripts.insert(assembled_transcript);
        self.transcript_models.extend(transcript_models);
    }

    pub fn get_size(&self) -> usize {
        self.assembled_transcripts.len()
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<Box<str>, u16>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<Box<str>, usize>) {
        self.read_names_map = read_names_map;
    }
}
