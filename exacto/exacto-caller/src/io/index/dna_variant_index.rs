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


use std::collections::{BTreeMap, HashMap, HashSet};

use crate::prelude::*;


pub struct DNAVariantIndex<'a> {
    records: &'a [DNAVariantRecord],

    /// HashMap<dna_variant_id, records index>
    by_variant_id: HashMap<u32, usize>,

    /// HashMap<chromosome_1, BTreeMap<position_1, Vec<records index>>>
    by_chromosome_1: HashMap<Box<str>, BTreeMap<u32, Vec<usize>>>,

    /// HashMap<chromosome_2, BTreeMap<position_2, Vec<records index>>>
    by_chromosome_2: HashMap<Box<str>, BTreeMap<u32, Vec<usize>>>
}

impl<'a> DNAVariantIndex<'a> {
    pub fn new(records: &'a [DNAVariantRecord]) -> Self {
        let mut by_variant_id: HashMap<u32, usize> = HashMap::new();
        let mut by_chromosome_1: HashMap<Box<str>, BTreeMap<u32, Vec<usize>>> = HashMap::new();
        let mut by_chromosome_2: HashMap<Box<str>, BTreeMap<u32, Vec<usize>>> = HashMap::new();

        for (i, r) in records.iter().enumerate() {
            if by_variant_id.insert(r.variant_id, i).is_some() {
                panic!(
                    "duplicate variant_id {} in records — variant_id must be unique",
                    r.variant_id
                );
            }

            by_chromosome_1
                .entry(r.chromosome_1.clone())
                .or_default()
                .entry(r.position_1)
                .or_default()
                .push(i);

            by_chromosome_2
                .entry(r.chromosome_2.clone())
                .or_default()
                .entry(r.position_2)
                .or_default()
                .push(i);
        }

        Self {
            records,
            by_variant_id,
            by_chromosome_1,
            by_chromosome_2
        }
    }

    pub fn records(&self) -> &'a [DNAVariantRecord] {
        self.records
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn get_by_variant_id(&self, variant_id: u32) -> Option<&'a DNAVariantRecord> {
        self.by_variant_id.get(&variant_id).map(|&i| &self.records[i])
    }

    pub fn get_by_range(
        &self,
        chromosome: &str,
        lo: u32,
        hi: u32
    ) -> Vec<&'a DNAVariantRecord> {
        let mut seen: HashSet<usize> = HashSet::new();
        let mut out: Vec<&'a DNAVariantRecord> = Vec::new();
        for map in [&self.by_chromosome_1, &self.by_chromosome_2] {
            if let Some(btree) = map.get(chromosome) {
                for (_pos, indices) in btree.range(lo..=hi) {
                    for &i in indices {
                        if seen.insert(i) {
                            out.push(&self.records[i]);
                        }
                    }
                }
            }
        }
        out
    }
}
