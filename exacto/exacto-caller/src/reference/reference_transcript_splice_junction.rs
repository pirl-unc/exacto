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
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};


#[derive(Debug,Serialize,Deserialize)]
pub struct ReferenceTranscriptSpliceJunction {
    pub chromosome: ReferenceChromosomeName,
    pub start: ReferencePosition,
    pub end: ReferencePosition,
    pub strand: Strand
}

impl PartialEq for ReferenceTranscriptSpliceJunction {
    fn eq(&self, other: &Self) -> bool {
        self.chromosome == other.chromosome &&
            self.start == other.start &&
            self.end == other.end &&
            self.strand == other.strand
    }
}

impl Eq for ReferenceTranscriptSpliceJunction {}

impl Hash for ReferenceTranscriptSpliceJunction {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.strand.hash(state);
    }
}

impl ReferenceTranscriptSpliceJunction {
    pub fn new(
        chromosome: &str,
        start: ReferencePosition,
        end: ReferencePosition,
        strand: Strand
    ) -> Self {
        assert!(
            start < end,
            "start should be less than end: {}-{}",
            start,
            end
        );

        Self {
            chromosome: chromosome.into(),
            start: start,
            end: end,
            strand: strand
        }
    }
}

impl Clone for ReferenceTranscriptSpliceJunction {
    fn clone(&self) -> Self {
        ReferenceTranscriptSpliceJunction {
            chromosome: self.chromosome.clone(),
            start: self.start,
            end: self.end,
            strand: self.strand.clone()
        }
    }
}
