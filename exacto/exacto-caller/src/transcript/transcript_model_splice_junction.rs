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
use std::hash::{Hash, Hasher};
use serde::{Deserialize, Serialize};


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModelSpliceJunction {
    pub reference_chromosome_id_1: ReferenceChromosomeID,
    pub reference_chromosome_id_2: ReferenceChromosomeID,
    pub reference_position_1: ReferencePosition,
    pub reference_position_2: ReferencePosition,
    pub reference_strand_1: Strand,
    pub reference_strand_2: Strand,
    pub splice_junction_number: SpliceJunctionNumber,
    pub read_position_1: ReadPosition,           // FASTX read sequence position immediately before the splice junction starts
    pub read_position_2: ReadPosition            // FASTX read sequence position immediately after the splice junction ends
}

impl PartialEq for TranscriptModelSpliceJunction {
    fn eq(&self, other: &Self) -> bool {
        self.reference_chromosome_id_1 == other.reference_chromosome_id_1 &&
            self.reference_chromosome_id_2 == other.reference_chromosome_id_2 &&
            self.reference_position_1 == other.reference_position_1 &&
            self.reference_position_2 == other.reference_position_2 &&
            self.splice_junction_number == other.splice_junction_number &&
            self.reference_strand_1 == other.reference_strand_1 &&
            self.reference_strand_2 == other.reference_strand_2 &&
            self.read_position_1 == other.read_position_1 &&
            self.read_position_2 == other.read_position_2
    }
}

impl Eq for TranscriptModelSpliceJunction {}

impl Hash for TranscriptModelSpliceJunction {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_chromosome_id_1.hash(state);
        self.reference_chromosome_id_2.hash(state);
        self.reference_position_1.hash(state);
        self.reference_position_2.hash(state);
        self.splice_junction_number.hash(state);
        self.reference_strand_1.hash(state);
        self.reference_strand_2.hash(state);
        self.read_position_1.hash(state);
        self.read_position_2.hash(state);
    }
}

impl TranscriptModelSpliceJunction {
    pub fn new(
        reference_chromosome_id_1: ReferenceChromosomeID,
        reference_chromosome_id_2: ReferenceChromosomeID,
        reference_position_1: ReferencePosition,
        reference_position_2: ReferencePosition,
        reference_strand_1: Strand,
        reference_strand_2: Strand,
        splice_junction_number: SpliceJunctionNumber,
        read_position_1: ReadPosition,
        read_position_2: ReadPosition
    ) -> Self {
        assert!(
            read_position_1 < read_position_2,
            "read_position_2 should be a read position following read_position_1: {}-{}",
            read_position_1,
            read_position_2
        );

        Self {
            reference_chromosome_id_1: reference_chromosome_id_1,
            reference_chromosome_id_2: reference_chromosome_id_2,
            reference_position_1: reference_position_1,
            reference_position_2: reference_position_2,
            reference_strand_1: reference_strand_1.clone(),
            reference_strand_2: reference_strand_2.clone(),
            splice_junction_number: splice_junction_number,
            read_position_1: read_position_1,
            read_position_2: read_position_2
        }
    }
}

impl Clone for TranscriptModelSpliceJunction {
    fn clone(&self) -> Self {
        TranscriptModelSpliceJunction {
            reference_chromosome_id_1: self.reference_chromosome_id_1,
            reference_chromosome_id_2: self.reference_chromosome_id_2,
            reference_position_1: self.reference_position_1,
            reference_position_2: self.reference_position_2,
            reference_strand_1: self.reference_strand_1.clone(),
            reference_strand_2: self.reference_strand_2.clone(),
            splice_junction_number: self.splice_junction_number,
            read_position_1: self.read_position_1,
            read_position_2: self.read_position_2
        }
    }
}
