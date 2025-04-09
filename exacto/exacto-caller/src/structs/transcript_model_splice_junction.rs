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
pub struct TranscriptModelSpliceJunction {
    pub chromosome_id: u16,
    pub start: u32,
    pub end: u32,
    pub number: u16,
    pub splice_signal_start: Box<str>,
    pub splice_signal_end: Box<str>,
    pub strand: Strands
}

impl PartialEq for TranscriptModelSpliceJunction {
    fn eq(&self, other: &Self) -> bool {
        self.chromosome_id == other.chromosome_id &&
            self.start == other.start &&
            self.end == other.end &&
            self.number == other.number &&
            self.splice_signal_start == other.splice_signal_start &&
            self.splice_signal_end == other.splice_signal_end &&
            self.strand == other.strand
    }
}

impl Eq for TranscriptModelSpliceJunction {}

impl Hash for TranscriptModelSpliceJunction {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome_id.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.number.hash(state);
        self.splice_signal_start.hash(state);
        self.splice_signal_end.hash(state);
        self.strand.hash(state);
    }
}

impl TranscriptModelSpliceJunction {
    pub fn new(
        chromosome_id: u16,
        start: u32,
        end: u32,
        number: u16,
        splice_signal_start: &str,
        splice_signal_end: &str,
        strand: Strands
    ) -> Self {
        Self {
            chromosome_id: chromosome_id,
            start: start,
            end: end,
            number: number,
            splice_signal_start: splice_signal_start.into(),
            splice_signal_end: splice_signal_end.into(),
            strand: strand.clone()
        }
    }
}

impl Clone for TranscriptModelSpliceJunction {
    fn clone(&self) -> Self {
        TranscriptModelSpliceJunction {
            chromosome_id: self.chromosome_id,
            start: self.start,
            end: self.end,
            number: self.number,
            splice_signal_start: self.splice_signal_start.clone(),
            splice_signal_end: self.splice_signal_end.clone(),
            strand: self.strand.clone()
        }
    }
}
