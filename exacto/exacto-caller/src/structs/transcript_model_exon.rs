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


use serde::{Deserialize, Serialize};
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};

use crate::common::constants::Strands;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModelExon {
    pub chromosome_id: u16,
    pub start: u32,
    pub end: u32,
    pub number: u16,
    pub strand: Strands
}

impl PartialEq for TranscriptModelExon {
    fn eq(&self, other: &Self) -> bool {
        self.chromosome_id == other.chromosome_id &&
            self.start == other.start &&
            self.end == other.end &&
            self.number == other.number &&
            self.strand == other.strand
    }
}

impl Eq for TranscriptModelExon {}

impl Hash for TranscriptModelExon {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome_id.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.number.hash(state);
        self.strand.hash(state);
    }
}

impl TranscriptModelExon {
    pub fn new(
        chromosome_id: u16,
        start: u32,
        end: u32,
        number: u16,
        strand: Strands
    ) -> Self {
        Self {
            chromosome_id: chromosome_id,
            start: start,
            end: end,
            number: number,
            strand: strand.clone()
        }
    }
}

impl Clone for TranscriptModelExon {
    fn clone(&self) -> Self {
        TranscriptModelExon {
            chromosome_id: self.chromosome_id,
            start: self.start,
            end: self.end,
            number: self.number,
            strand: self.strand.clone()
        }
    }
}
