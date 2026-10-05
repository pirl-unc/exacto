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
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::prelude::{SpliceJunction, SplicingEventReadCounts};


pub struct RNANovelJunctionReadCounts {
    /// Vec<(chromosome ID, first intron base, last intron base)>, sorted.
    novel_junctions: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)>,

    num_spliced: Vec<ReadSupport>,
    depths: Vec<ReadDepth>
}

impl RNANovelJunctionReadCounts {
    pub fn new(
        splice_junctions: &[Vec<SpliceJunction>],
        is_annotated: impl Fn(&(ReferenceChromosomeID, ReferencePosition, ReferencePosition)) -> bool
    ) -> Self {
        // Step 1. Count the reads splicing each intron.
        // HashMap<(chromosome ID, first intron base, last intron base), number of reads>
        let mut spliced_counts: HashMap<(ReferenceChromosomeID, ReferencePosition, ReferencePosition), ReadSupport> = HashMap::new();
        for read_splice_junctions in splice_junctions.iter() {
            for intron in identify_spliced_introns(read_splice_junctions) {
                *spliced_counts.entry(intron).or_insert(0) += 1;
            }
        }
        let mut introns: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = spliced_counts.keys().copied().collect();
        introns.sort_unstable();
        let is_novel: Vec<bool> = introns.iter().map(|intron| !is_annotated(intron)).collect();

        // Step 2. The introns are sorted by their first base, so the introns overlapping one are
        // found walking down from the last intron that starts at or before its last base, until
        // no intron further down reaches its first base.
        // Vec<last base of the introns of the chromosome up to this one>
        let mut max_intron_ends: Vec<ReferencePosition> = Vec::with_capacity(introns.len());
        for (index, &(chromosome, _, intron_end)) in introns.iter().enumerate() {
            if index > 0 && introns[index - 1].0 == chromosome {
                max_intron_ends.push(max_intron_ends[index - 1].max(intron_end));
            } else {
                max_intron_ends.push(intron_end);
            }
        }
        let identify_overlapping_introns = |(chromosome, intron_start, intron_end): (ReferenceChromosomeID, ReferencePosition, ReferencePosition)| -> Vec<usize> {
            let last: usize = introns.partition_point(|&other| other <= (chromosome, intron_end, ReferencePosition::MAX));
            let mut overlapping: Vec<usize> = Vec::new();
            for index in (0..last).rev() {
                if introns[index].0 != chromosome || max_intron_ends[index] < intron_start {
                    break;
                }
                if introns[index].2 >= intron_start {
                    overlapping.push(index);
                }
            }
            overlapping
        };

        // Step 3. Count the reads splicing each novel splice junction or an intron overlapping it.
        let mut depths: Vec<ReadDepth> = vec![0; introns.len()];
        for read_splice_junctions in splice_junctions.iter() {
            let overlapped: BTreeSet<usize> = identify_spliced_introns(read_splice_junctions)
                .into_iter()
                .flat_map(|intron| identify_overlapping_introns(intron))
                .filter(|&index| is_novel[index])
                .collect();
            for index in overlapped {
                depths[index] += 1;
            }
        }

        // Step 4. Keep the novel splice junctions.
        let novel: Vec<usize> = (0..introns.len()).filter(|&index| is_novel[index]).collect();
        Self {
            novel_junctions: novel.iter().map(|&index| introns[index]).collect(),
            num_spliced: novel.iter().map(|&index| spliced_counts[&introns[index]]).collect(),
            depths: novel.iter().map(|&index| depths[index]).collect()
        }
    }

    pub fn get_depths(&self) -> HashSet<ReadDepth> {
        self.depths.iter().copied().collect()
    }
}

impl SplicingEventReadCounts for RNANovelJunctionReadCounts {
    type Input = Vec<SpliceJunction>;
    
    fn get_read_counts(&self, splice_junctions: &Vec<SpliceJunction>) -> Vec<(ReadSupport, ReadDepth)> {
        identify_spliced_introns(splice_junctions)
            .into_iter()
            .filter_map(|intron| self.novel_junctions.binary_search(&intron).ok())
            .map(|index| (self.num_spliced[index], self.depths[index]))
            .collect()
    }
}


fn identify_spliced_introns(splice_junctions: &[SpliceJunction]) -> BTreeSet<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> {
    splice_junctions
        .iter()
        .filter(|junction| junction.chromosome_1 == junction.chromosome_2)
        .map(|junction| junction.intron_span_key())
        .collect()
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/novel_junction_read_counts.rs"]
mod tests;