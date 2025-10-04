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
            variant_records: HashSet::new()
        }
    }

    pub fn add_variant_record(&mut self, variant_record: VariantRecord) {
        self.variant_records.insert(variant_record);
    }
    
    /// Get the consensus VariantRecord object for this VariantCall object.
    ///
    /// # Returns
    ///
    /// (VariantRecord,read IDs).
    pub fn get_consensus_record(&self) -> (&VariantRecord, Vec<usize>) {
        let mut map: HashMap<(u16, u32, GraphOperationType, u16, u32, GraphOperationType, Box<str>, VariantType),Vec<&VariantRecord>> = HashMap::new();
        for variant_record in self.variant_records.iter() {
            let key: (u16, u32, GraphOperationType, u16, u32, GraphOperationType, Box<str>, VariantType) = (
                variant_record.graph_operation.get_chromosome_1(),
                variant_record.graph_operation.get_position_1(),
                variant_record.graph_operation.get_operation_type_1().clone(),
                variant_record.graph_operation.get_chromosome_2(),
                variant_record.graph_operation.get_position_2(),
                variant_record.graph_operation.get_operation_type_2().clone(),
                variant_record.get_standardized_sequence().into(),
                variant_record.graph_operation.get_variant_type().clone()
            );
            map
                .entry(key)
                .or_insert(Vec::new())
                .push(variant_record);
        }
        let max_vec = map
            .iter()
            .max_by_key(|(_, v)| v.len())
            .expect("self.variant_records is empty.")
            .1;
        (max_vec[0], max_vec.iter().map(|v| v.read_id).collect())
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
            variant_records: self.variant_records.clone()
        }
    }
}
