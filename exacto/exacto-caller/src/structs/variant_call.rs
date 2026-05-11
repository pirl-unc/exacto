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
use exacto_core::prelude::{reverse_complement, Strand};
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Eq,Serialize,Deserialize)]
pub struct VariantCall {
    pub id: usize,
    pub total_depth: i32,
    pub variant_records: HashSet<VariantRecord>
}

impl Hash for VariantCall {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let variant_record: &VariantRecord = self.get_consensus_record().0;
        variant_record.graph_operation.get_chromosome_1().hash(state);
        variant_record.graph_operation.get_position_1().hash(state);
        variant_record.graph_operation.get_operation_type_1().hash(state);
        variant_record.graph_operation.get_chromosome_2().hash(state);
        variant_record.graph_operation.get_position_2().hash(state);
        variant_record.graph_operation.get_operation_type_2().hash(state);
        variant_record.graph_operation.get_variant_type().hash(state);
        let mut read_ids: Vec<usize> = Vec::new();
        for variant_record in self.variant_records.iter() {
            read_ids.push(variant_record.read_id);
        }
        read_ids.sort();
        for read_id in read_ids.iter() {
            read_id.hash(state);
        }
    }
}

impl PartialEq for VariantCall {
    fn eq(&self, other: &Self) -> bool {
        let variant_record_1: &VariantRecord = self.get_consensus_record().0;
        let variant_record_2: &VariantRecord = other.get_consensus_record().0;
        if variant_record_1.get_chromosome_1() == variant_record_2.get_chromosome_1() &&
            variant_record_1.get_position_1() == variant_record_2.get_position_1() &&
            variant_record_1.get_operation_1() == variant_record_2.get_operation_1() &&
            variant_record_1.get_chromosome_2() == variant_record_2.get_chromosome_2() &&
            variant_record_1.get_position_2() == variant_record_2.get_position_2() &&
            variant_record_1.get_operation_2() == variant_record_2.get_operation_2() &&
            variant_record_1.get_variant_type() == variant_record_2.get_variant_type() {
            let mut read_ids_1: Vec<usize> = self.get_read_ids();
            let mut read_ids_2: Vec<usize> = other.get_read_ids();
            read_ids_1.sort();
            read_ids_2.sort();
            if read_ids_1 == read_ids_2 {
                true
            } else {
                false
            }
        } else {
            false
        }
    }
}

impl VariantCall {
    pub fn new(id: usize) -> Self {
        Self {
            id: id,
            total_depth: -1,
            variant_records: HashSet::new()
        }
    }

    pub fn add_variant_record(&mut self, variant_record: VariantRecord) {
        self.variant_records.insert(variant_record);
    }

    pub fn get_alternate_allele_fraction(&self) -> f32 {
        assert!(self.total_depth > 0);
        let num_alt_allele_read_count: usize = self.get_read_ids().len();
        num_alt_allele_read_count as f32 / self.total_depth as f32
    }
    
    /// Get the consensus VariantRecord object for this VariantCall object.
    ///
    /// # Returns
    ///
    /// (VariantRecord,read IDs).
    pub fn get_consensus_record(&self) -> (&VariantRecord, Vec<usize>) {
        assert!(!self.variant_records.is_empty(), "self.variant_records is empty");

        type Key = (
            u16,
            u32,
            GraphOperationType,
            u16,
            u32,
            GraphOperationType,
            Box<str>,
            VariantType
        );

        // Build a stable "key" function
        let make_key = |vr: &VariantRecord| -> Key {
            (
                vr.graph_operation.get_chromosome_1(),
                vr.graph_operation.get_position_1(),
                vr.graph_operation.get_operation_type_1().clone(),
                vr.graph_operation.get_chromosome_2(),
                vr.graph_operation.get_position_2(),
                vr.graph_operation.get_operation_type_2().clone(),
                vr.get_standardized_sequence().into(),
                vr.graph_operation.get_variant_type().clone()
            )
        };

        let mut map: HashMap<Key, Vec<&VariantRecord>> = HashMap::new();

        for vr in self.variant_records.iter() {
            let key: Key = make_key(vr);
            map.entry(key).or_insert_with(Vec::new).push(vr);
        }

        // If there is a consensus group (>=2), return it
        if let Some((_key, max_vec)) = map.iter().max_by_key(|(_, v)| v.len()) {
            if max_vec.len() >= 2 {
                let consensus_record = max_vec[0];
                let read_ids = max_vec.iter().map(|v| v.read_id).collect();
                return (consensus_record, read_ids);
            }
        }

        // Otherwise, choose the record with the median variant size.
        // Then tie-break by highest support among records with that same size.
        let mut sizes: Vec<isize> = self.variant_records.iter().map(|vr| vr.get_variant_size()).collect();
        sizes.sort_unstable();
        let median_size: isize = sizes[(sizes.len() - 1) / 2];

        // Filter candidates that match median size
        let mut median_candidates: Vec<&VariantRecord> = self
            .variant_records
            .iter()
            .filter(|vr| vr.get_variant_size() == median_size)
            .collect();
        assert!(!median_candidates.is_empty(), "median_candidates is empty");

        // For insertion edge-case: same size, different sequences.
        // Pick the candidate whose key-group has the largest support count.
        // Tie-break deterministically by smallest read_id.
        median_candidates.sort_by(|a, b| {
            let ka: Key = make_key(a);
            let kb: Key = make_key(b);

            let sa: usize = map.get(&ka).map(|v| v.len()).unwrap_or(0);
            let sb: usize = map.get(&kb).map(|v| v.len()).unwrap_or(0);

            // Primary: larger support first
            sb.cmp(&sa)
                // Tie-breaker: smaller read_id first (deterministic)
                .then_with(|| a.read_id.cmp(&b.read_id))
        });

        let chosen: &VariantRecord = median_candidates[0];
        let chosen_key: Key = make_key(chosen);

        let read_ids: Vec<usize> = map
            .get(&chosen_key)
            .map(|v| v.iter().map(|vr| vr.read_id).collect())
            .unwrap_or_else(|| vec![chosen.read_id]);

        (chosen, read_ids)
    }

    pub fn get_named_consensus_record(&self, read_names_map: &BiMap<Box<str>,usize>) -> (&VariantRecord, Vec<Box<str>>) {
        let (consensus_variant_record, read_ids) = self.get_consensus_record();
        let mut read_names: Vec<Box<str>> = Vec::new();
        for read_id in read_ids.iter() {
            read_names.push(read_names_map.get_by_right(&read_id).unwrap().clone());
        }
        (consensus_variant_record, read_names)
    }

    pub fn get_read_ids(&self) -> Vec<usize> {
        self.variant_records
            .iter()
            .map(|record| record.read_id.clone())
            .collect()
    }

    pub fn get_read_names(&self, read_names_map: &HashMap<usize,Box<str>>) -> Vec<Box<str>> {
        self.variant_records
            .iter()
            .map(|record| read_names_map.get(&record.read_id).unwrap().clone())
            .collect()
    }
    
    pub fn get_sequence_operation_boxed_str(&self) -> Vec<Box<str>> {
        self.variant_records
            .iter()
            .map(|record| record.get_graph_operation_boxed_str())
            .collect()
    }

    pub fn get_sequence_operation_named_boxed_str(&self, chromosome_names_map: &BiMap<Box<str>,u16>) -> Vec<Box<str>> {
        self.variant_records
            .iter()
            .map(|record| record.get_graph_operation_named_boxed_str(chromosome_names_map))
            .collect()
    }

    pub fn get_strand_1(&self) -> Strand {
        let has_forward: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.graph_operation.get_strand_1() == &Strand::Forward);
        let has_reverse: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.graph_operation.get_strand_1() == &Strand::Reverse);
        match (has_forward, has_reverse) {
            (true, true)    => Strand::Both,
            (true, false)   => Strand::Forward,
            (false, true)   => Strand::Reverse,
            (false, false)  => panic!("No strand information available for this variant call.")
        }
    }

    pub fn get_strand_2(&self) -> Strand {
        let has_forward: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.graph_operation.get_strand_2() == &Strand::Forward);
        let has_reverse: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.graph_operation.get_strand_2() == &Strand::Reverse);
        match (has_forward, has_reverse) {
            (true, true)    => Strand::Both,
            (true, false)   => Strand::Forward,
            (false, true)   => Strand::Reverse,
            (false, false)  => panic!("No strand information available for this variant call.")
        }
    }
    
    pub fn get_total_depth(&self) -> i32 {
        self.total_depth
    }

    pub fn set_total_depth(&mut self, total_depth: i32) {
        self.total_depth = total_depth;
    }

    pub fn to_tsv_string(
        &self,
        chromosome_names_map: &BiMap<Box<str>,u16>,
        read_names_map: &BiMap<Box<str>,usize>
    ) -> String {
        let (consensus_record, consensus_read_names) = self.get_named_consensus_record(read_names_map);
        let chromosome_1: &str = &*chromosome_names_map.get_by_right(&consensus_record.get_chromosome_1()).unwrap();
        let chromosome_2: &str = &*chromosome_names_map.get_by_right(&consensus_record.get_chromosome_2()).unwrap();
        let read_names: Vec<&str> = self.get_read_ids()
            .iter()
            .map(|read_id| &**read_names_map.get_by_right(read_id).unwrap())
            .collect();
        let sequence: String = if consensus_record.get_strand_1().clone() == Strand::Forward {
            consensus_record.graph_operation.get_sequence().to_string().to_uppercase()
        } else {
            reverse_complement(&*consensus_record.graph_operation.get_sequence()).to_string().to_uppercase()
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.id,
            chromosome_1,
            consensus_record.graph_operation.get_position_1(),
            consensus_record.graph_operation.get_strand_1().as_str(),
            consensus_record.graph_operation.get_operation_type_1().as_str(),
            chromosome_2,
            consensus_record.graph_operation.get_position_2(),
            consensus_record.graph_operation.get_strand_2().as_str(),
            consensus_record.graph_operation.get_operation_type_2().as_str(),
            consensus_record.get_variant_size() as i64,
            consensus_record.get_variant_type().as_str().to_string(),
            sequence,
            consensus_read_names.join(","),
            consensus_read_names.len() as u64,
            read_names.join(","),
            read_names.len() as u64
        )
    }
}

impl Clone for VariantCall {
    fn clone(&self) -> Self {
        VariantCall {
            id: self.id,
            total_depth: self.total_depth,
            variant_records: self.variant_records.clone()
        }
    }
}
