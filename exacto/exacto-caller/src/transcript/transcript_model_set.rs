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
use exacto_core::prelude::{ReadID, ReadName, ReferenceChromosomeID, ReferenceChromosomeName};
use serde::{Deserialize, Serialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptModelSet {
    pub transcript_models: Vec<TranscriptModel>,
    pub read_names_map: BiMap<ReadName, ReadID>,
    pub chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID>
}

impl TranscriptModelSet {
    pub fn new() -> Self {
        Self {
            transcript_models: Vec::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add(&mut self, transcript_model: TranscriptModel) {
        self.transcript_models.push(transcript_model);
    }

    pub fn get_size(&self) -> usize {
        self.transcript_models.len()
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<ReadName, ReadID>) {
        self.read_names_map = read_names_map;
    }
}