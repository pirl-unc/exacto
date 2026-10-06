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


use exacto_core::prelude::{ReadID, ReferenceChromosomeID, ReferencePosition};
use serde::{Deserialize, Serialize};

use crate::prelude::SpliceJunction;


#[derive(Debug, Serialize, Deserialize)]
pub struct RNAReadCharacterizationSummary {
    pub read_id: ReadID,

    /// Chromosome of the model's reference span; `None` when start and end sit on different
    /// chromosomes, which the clusterer drops anyway.
    pub chromosome: Option<ReferenceChromosomeID>,

    /// Reference positions of the model's first and last exon, in read order.
    pub reference_start: ReferencePosition,
    pub reference_end: ReferencePosition,

    /// Whether the read ends at `reference_start` and at `reference_end`. It does not where
    /// read bases the aligner did not place lie beyond the position.
    pub ends_at_reference_start: bool,
    pub ends_at_reference_end: bool,

    pub splice_junctions: Vec<SpliceJunction>,
    pub exons: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)>,

    /// Start of the serialized model.
    pub offset: u64,

    /// Serialized byte length.
    pub length: u64
}


impl RNAReadCharacterizationSummary {
    pub fn identify_retained_introns(
        &self,
        introns: &[(ReferenceChromosomeID, ReferencePosition, ReferencePosition)],
        introns_by_end: &[(ReferenceChromosomeID, ReferencePosition, ReferencePosition)]
    ) -> Vec<usize> {
        let mut retained: Vec<usize> = Vec::new();
        for &(chromosome, exon_start, exon_end) in self.exons.iter() {
            let ends_at_exon_start: bool = self.ends_at_reference_start
                && self.chromosome == Some(chromosome)
                && self.reference_start == exon_start;
            let ends_at_exon_end: bool = self.ends_at_reference_end
                && self.chromosome == Some(chromosome)
                && self.reference_end == exon_end;

            // The introns starting inside the block. The block goes on past the intron, or the
            // read ends here.
            let first: usize = introns.partition_point(|&intron| intron <= (chromosome, exon_start, u32::MAX));
            for index in first..introns.len() {
                let (intron_chromosome, intron_start, intron_end): (ReferenceChromosomeID, ReferencePosition, ReferencePosition) = introns[index];
                if intron_chromosome != chromosome || intron_start > exon_end {
                    break;
                }
                if intron_end < exon_end || ends_at_exon_end {
                    retained.push(index);
                }
            }

            // The introns ending inside the block that the read starts among the bases of.
            if !ends_at_exon_start {
                continue;
            }
            let first: usize = introns_by_end.partition_point(|&intron| intron < (chromosome, exon_start, 0));
            for &(intron_chromosome, intron_end, intron_start) in introns_by_end[first..].iter() {
                if intron_chromosome != chromosome || intron_end >= exon_end {
                    break;
                }
                if intron_start > exon_start {
                    continue;
                }
                if let Ok(index) = introns.binary_search(&(chromosome, intron_start, intron_end)) {
                    retained.push(index);
                }
            }
        }

        // A read that splices an intron in one record and covers it in another spliced it.
        retained.retain(|&index| {
            self.splice_junctions
                .iter()
                .all(|junction| junction.intron_span_key() != introns[index])
        });
        retained.sort_unstable();
        retained.dedup();
        retained
    }
}
