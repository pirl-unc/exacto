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


use exacto_core::prelude::{ReferenceChromosomeName, ReferencePosition};
use std::collections::{BTreeMap, HashMap};


pub(crate) fn generate_genomic_regions(
    regions: &Vec<(&str, ReferencePosition, ReferencePosition)>,
    chromosomes: &Vec<ReferenceChromosomeName>,
    chromosome_lengths: &HashMap<ReferenceChromosomeName, u32>,
    chunk_size: u32
) -> BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> {
    let ordered_regions: BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = if regions.is_empty() {
        let regions: HashMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = exacto_core::prelude::generate_regions(
            &chromosomes,
            &chromosome_lengths,
            chunk_size
        );
        let mut ordered_regions: BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = regions.into_iter().collect();
        for (_, vec) in ordered_regions.iter_mut() {
            vec.sort_by(|a, b| a.0.cmp(&b.0));
        }
        ordered_regions
    } else {
        let regions_split: Vec<(ReferenceChromosomeName, ReferencePosition, ReferencePosition)> = exacto_core::prelude::split_regions(regions, chunk_size);
        let mut ordered_regions: BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = BTreeMap::new();
        for (contig, start, end) in regions_split {
            ordered_regions
                .entry(contig)
                .or_insert_with(Vec::new)
                .push((start, end));
        }
        for vec in ordered_regions.values_mut() {
            vec.sort_by_key(|(start, _)| *start);
        }
        ordered_regions
    };
    
    ordered_regions
}


pub(crate) fn get_region_ends(chunks: &[(ReferencePosition, ReferencePosition)]) -> Vec<ReferencePosition> {
    let mut region_ends: Vec<ReferencePosition> = chunks.iter().map(|&(_, end)| end).collect();
    for i in (1..chunks.len()).rev() {
        if chunks[i].0 <= chunks[i - 1].1.saturating_add(1) {
            region_ends[i - 1] = region_ends[i - 1].max(region_ends[i]);
        }
    }
    region_ends
}
