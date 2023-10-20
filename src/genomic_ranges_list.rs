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


use crate::genomic_range::GenomicRange;
use std::collections::HashMap;


#[derive(Clone, Debug)]
pub struct GenomicRangesList {
    genomic_ranges: HashMap<String, Vec<GenomicRange>> // key = chromosome
}

impl GenomicRangesList {
    fn new() -> Self {
        GenomicRangesList {
            genomic_ranges: HashMap::new()
        }
    }

    pub fn add_genomic_range(&mut self, genomic_range: GenomicRange) {
        let key = genomic_range.chromosome.clone();

        // Find an index position to insert the new GenomicRange object
        let insert_idx = self.genomic_ranges.entry(key.clone()).or_insert(Vec::new()).binary_search_by(|item| item.start.cmp(&genomic_range.start));
        match insert_idx {
            Ok(idx) => {
                self.genomic_ranges.entry(key).or_insert(Vec::new()).insert(idx, genomic_range);
            }
            Err(idx) => {
                self.genomic_ranges.entry(key).or_insert(Vec::new()).insert(idx, genomic_range);
            }
        }
    }

    pub fn find_overlaps(&mut self, chromosome: String, start: isize, end: isize) -> Vec<GenomicRange> {
        let mut overlaps: Vec<GenomicRange> = Vec::new();
        let key = chromosome.clone();
        for genomic_ranges in self.genomic_ranges.get(&key) {
            for genomic_range in genomic_ranges {
                let mut genomic_range_ = genomic_range.clone();
                if genomic_range_.overlaps(chromosome.clone(), start, end) {
                    overlaps.push(genomic_range_);
                }
            }
        }
        return overlaps;
    }
}
