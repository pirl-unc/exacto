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


use bimap::BiMap;
use exacto_core::prelude::*;
use polars::prelude::*;
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::str::FromStr;
use rayon::ThreadPool;
use tempfile::TempPath;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct DNAVariantCallSet {
    pub origin: DNAVariantOrigin,
    
    pub variant_calls: HashMap<VariantID, VariantCall>,

    /// BiMap<read name, read ID>
    pub read_names_map: BiMap<ReadName, ReadID>,

    /// BiMap<chromosome name, chromosome ID>
    pub chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,

    /// HashMap<chromsoome ID, BTreeMap<position, HashSet<variant call ID>>>
    position_index: HashMap<ReferenceChromosomeID, BTreeMap<ReferencePosition, HashSet<VariantID>>>
}

impl DNAVariantCallSet {
    pub fn new(origin: DNAVariantOrigin) -> Self {
        Self {
            origin: origin,
            variant_calls: HashMap::new(),
            position_index: HashMap::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add_variant_call(&mut self, variant_call: VariantCall) {
        if self.variant_calls.contains_key(&variant_call.get_id()) == false {
            let consensus_graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
            let chromosome_1_id: ReferenceChromosomeID = consensus_graph_operation.get_chromosome_1();
            let chromosome_2_id: ReferenceChromosomeID = consensus_graph_operation.get_chromosome_2();
            let position_1: ReferencePosition = consensus_graph_operation.get_position_1();
            let position_2: ReferencePosition = consensus_graph_operation.get_position_2();

            // Index the variant by its first breakpoint position.
            self.position_index
                .entry(chromosome_1_id)
                .or_insert_with(BTreeMap::new)
                .entry(position_1)
                .or_insert_with(HashSet::new)
                .insert(variant_call.get_id());

            // Index the variant by its second breakpoint position.
            self.position_index
                .entry(chromosome_2_id)
                .or_insert_with(BTreeMap::new)
                .entry(position_2)
                .or_insert_with(HashSet::new)
                .insert(variant_call.get_id());

            // Add to self.variant_calls.
            self.variant_calls.insert(variant_call.get_id(), variant_call);
        }
    }

    pub fn get_size(&self) -> usize {
        self.variant_calls.len()
    }

    pub fn get_variant_calls(&self) -> Vec<&VariantCall> {
        let mut variant_calls: Vec<&VariantCall> = self.variant_calls.values().collect();
        variant_calls.sort_unstable_by_key(|variant_call| variant_call.get_id());
        variant_calls
    }

    pub fn get_variant_calls_by_range(&self, chromosome_id: ReferenceChromosomeID, start: ReferencePosition, end: ReferencePosition) -> Vec<&VariantCall> {
        let mut result_ids = HashSet::new();
        if let Some(position_map) = self.position_index.get(&chromosome_id) {
            for (_pos, variant_call_ids) in position_map.range(start..=end) {
                result_ids.extend(variant_call_ids.iter().cloned());
            }
        }
        result_ids.iter()
            .filter_map(|id| self.variant_calls.get(id))
            .collect()
    }

    pub fn get_variant_records(&self) -> Vec<&VariantRecord> {
        let mut variant_records: Vec<&VariantRecord> = Vec::new();
        for variant_call in self.variant_calls.values() {
            for variant_record in variant_call.get_variant_records().iter() {
                variant_records.push(variant_record);
            }
        }
        variant_records
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<ReferenceChromosomeName,ReferenceChromosomeID>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<ReadName,ReadID>) {
        self.read_names_map = read_names_map;
    }

    pub fn remove_variant_call(&mut self, variant_call: &VariantCall) {
        self.variant_calls.remove(&variant_call.get_id());
    }
}

impl Clone for DNAVariantCallSet {
    fn clone(&self) -> Self {
        DNAVariantCallSet {
            origin: self.origin.clone(),
            variant_calls: self.variant_calls.clone(),
            position_index: self.position_index.clone(),
            read_names_map: self.read_names_map.clone(),
            chromosome_names_map: self.chromosome_names_map.clone()
        }
    }
}

pub(crate) fn merge_temp_variant_call_sets(
    origin: DNAVariantOrigin,
    temp_files: &[TempPath],
    read_names_map: &BiMap<ReadName, ReadID>,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    num_threads: usize
) -> DNAVariantCallSet {
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let variant_calls: HashSet<VariantCall> = thread_pool.install(|| {
        temp_files
            .par_iter()
            .map(|temp_path| {
                let file = File::open(temp_path).unwrap();
                let mut reader = BufReader::new(file);
                let set: DNAVariantCallSet = bincode::deserialize_from(&mut reader).expect("Failed to deserialize data");
                let mut local_variant_calls: HashSet<VariantCall> = HashSet::with_capacity(set.variant_calls.len());
                for (_, vc) in set.variant_calls {
                    local_variant_calls.insert(vc);
                }
                local_variant_calls
            })
            .reduce(
                || HashSet::new(),
                |mut acc, local| {
                    // Reduce rehashing during merge
                    acc.reserve(local.len());
                    acc.extend(local);
                    acc
                }
            )
    });

    // Number the calls in genomic order, so that two runs of one input give the same ids.
    let mut variant_calls: Vec<VariantCall> = variant_calls.into_iter().collect();
    variant_calls.sort_by_cached_key(|variant_call| {
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        (
            operation.get_chromosome_1(),
            operation.get_position_1(),
            operation.get_chromosome_2(),
            operation.get_position_2(),
            operation.as_boxed_str(),
            variant_call.get_read_ids()
        )
    });

    let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new(origin);

    let mut variant_call_id: VariantID = 1;
    for mut variant_call in variant_calls {
        variant_call.set_id(variant_call_id);
        variant_call_set.add_variant_call(variant_call);
        variant_call_id += 1;
    }

    variant_call_set.load_read_names(read_names_map.clone());
    variant_call_set.load_chromosome_names(chromosome_names_map.clone());

    variant_call_set
}
