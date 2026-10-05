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


#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SpliceJunction {
    pub chromosome_1: ReferenceChromosomeID,
    pub chromosome_2: ReferenceChromosomeID,
    pub position_1: ReferencePosition,
    pub position_2: ReferencePosition,
    pub strand_1: Strand,
    pub strand_2: Strand
}

impl SpliceJunction {
    pub fn new(
        chromosome_1: ReferenceChromosomeID,
        chromosome_2: ReferenceChromosomeID,
        position_1: ReferencePosition,
        position_2: ReferencePosition,
        strand_1: Strand,
        strand_2: Strand
    ) -> Self {
        Self {
            chromosome_1,
            chromosome_2,
            position_1,
            position_2,
            strand_1,
            strand_2
        }
    }

    /// The intron as a genomic interval. Position order encodes the strand, so the
    /// intron is the low..=high of the two positions regardless of direction.
    pub fn intron_span(&self) -> (ReferencePosition, ReferencePosition) {
        (
            self.position_1.min(self.position_2),
            self.position_1.max(self.position_2)
        )
    }

    /// The intron keyed the way the annotation indexes it: (chromosome, low, high).
    pub fn intron_span_key(&self) -> (ReferenceChromosomeID, ReferencePosition, ReferencePosition) {
        let (start, end): (ReferencePosition, ReferencePosition) = self.intron_span();
        (self.chromosome_1, start, end)
    }
}