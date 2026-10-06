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


use abpoa_rs::AlignmentMode;
use bimap::BiMap;
use exacto_consensus::prelude::perform_partial_order_alignment;
use exacto_core::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Eq,Serialize,Deserialize)]
pub struct VariantCall {
    id: VariantID,
    total_depth: i32,
    variant_records: HashSet<VariantRecord>,
    consensus_graph_operation: GraphOperation,
    consensus_read_ids: HashSet<ReadID>,
    consensus_method: GraphOperationConsensusMethod
}

impl Hash for VariantCall {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.get_consensus_graph_operation().hash(state);
        self.get_consensus_method().hash(state);
        for read_id in self.get_read_ids().iter() {
            read_id.hash(state);
        }
    }
}

impl PartialEq for VariantCall {
    fn eq(&self, other: &Self) -> bool {
        if *self.get_consensus_graph_operation() == *other.get_consensus_graph_operation() {
            self.get_read_ids() == other.get_read_ids()
        } else {
            false
        }
    }
}

impl VariantCall {
    pub fn from_variant_records(
        id: VariantID,
        variant_records: HashSet<VariantRecord>,
        match_score: i32,
        mismatch_score: i32,
        gap_open_score: i32,
        gap_extend_score: i32
    ) -> Self {
        // Step 1. The graph operations and chromosomes should all agree.
        let mut pool: Vec<&VariantRecord> = variant_records
            .iter()
            .filter(|vr| vr.is_resolved())
            .collect();
        if pool.is_empty() {
            pool = variant_records.iter().collect();
        }

        // A HashSet has no order; sort so every downstream "first" means the lowest read id.
        pool.sort_by_key(|vr| (vr.get_read_id(), vr.get_read_position_1()));

        let graph_operation_type_1: GraphOperationType = pool.first().unwrap().get_operation_1().clone();
        let graph_operation_type_2: GraphOperationType = pool.first().unwrap().get_operation_2().clone();
        let chromosome_1: ReferenceChromosomeID = pool.first().unwrap().get_chromosome_1();
        let chromosome_2: ReferenceChromosomeID = pool.first().unwrap().get_chromosome_2();
        for variant_record in variant_records.iter() {
            let side_1: (ReferenceChromosomeID, &GraphOperationType) = (variant_record.get_chromosome_1(), variant_record.get_operation_1());
            if variant_record.is_resolved() {
                assert_eq!(side_1, (chromosome_1, &graph_operation_type_1));
                assert_eq!((variant_record.get_chromosome_2(), variant_record.get_operation_2()), (chromosome_2, &graph_operation_type_2));
            } else {
                // A clip is support at one breakend of the junction, which may be either side.
                assert!(side_1 == (chromosome_1, &graph_operation_type_1) || side_1 == (chromosome_2, &graph_operation_type_2));
            }
        }

        // Step 2. Determine the consensus graph operation.
        let (consensus_graph_operation, consensus_read_ids, consensus_method):
            (GraphOperation, HashSet<ReadID>, GraphOperationConsensusMethod) = Self::resolve_consensus_operation(
            &variant_records,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score
        );

        Self {
            id: id,
            total_depth: -1,
            variant_records: variant_records,
            consensus_graph_operation: consensus_graph_operation,
            consensus_read_ids: consensus_read_ids,
            consensus_method: consensus_method
        }
    }

    pub fn resolve_consensus_operation(
        variant_records: &HashSet<VariantRecord>,
        match_score: i32,
        mismatch_score: i32,
        gap_open_score: i32,
        gap_extend_score: i32
    ) -> (GraphOperation, HashSet<ReadID>, GraphOperationConsensusMethod) {
        // Step 1. Get resolved variant records if they exist.
        let resolved: Vec<&VariantRecord> = variant_records
            .iter()
            .filter(|vr| *vr.get_operation_2() != GraphOperationType::Noop)
            .collect();
        let mut variant_records_: Vec<&VariantRecord> = if resolved.is_empty() {
            variant_records.iter().collect()
        } else {
            resolved
        };

        // A HashSet has no order; sort so "first" means the lowest read id and the partial
        // order alignment takes its sequences in one order.
        variant_records_.sort_by_key(|vr| (
            vr.get_read_id(),
            vr.get_read_position_1(),
            vr.get_read_position_2(),
            vr.get_position_1(),
            vr.get_position_2()
        ));

        // Step 2. Return the one variant record if that is the only variant record in this variant call.
        if variant_records_.len() == 1 {
            return (
                variant_records_.first().unwrap().get_graph_operation().clone(),
                HashSet::from_iter(vec![variant_records_.first().unwrap().get_read_id()]),
                GraphOperationConsensusMethod::Modal
            );
        }

        // Step 3. Create a map of variant records by "keys":
        // (chromosome 1, position 1, GraphOperationType,
        //  chromosome 2, position 2, GraphOperationType,
        //  sequence, VariantType)
        type Key = (ReferenceChromosomeID, ReferencePosition, GraphOperationType, ReferenceChromosomeID, ReferencePosition, GraphOperationType, Box<str>, VariantType);

        // Build a stable "key" function
        let make_key = |vr: &VariantRecord| -> Key {
            (
                vr.get_chromosome_1(),
                vr.get_position_1(),
                vr.get_operation_1().clone(),
                vr.get_chromosome_2(),
                vr.get_position_2(),
                vr.get_operation_2().clone(),
                vr.get_standardized_sequence().into(),
                vr.get_variant_type().clone()
            )
        };

        let mut map: HashMap<Key, Vec<&VariantRecord>> = HashMap::new();

        for variant_record in variant_records_.iter() {
            let key: Key = make_key(variant_record);
            map
                .entry(key)
                .or_insert_with(Vec::new)
                .push(variant_record);
        }

        // Step 4. If there is a consensus group (>=2), return it.
        let mut groups: Vec<(&Key, &Vec<&VariantRecord>)> = map.iter().collect();
        groups.sort_by(|(key_a, vec_a), (key_b, vec_b)| {
            vec_b.len().cmp(&vec_a.len())
                .then_with(|| key_a.1.cmp(&key_b.1))   // position 1
                .then_with(|| key_a.4.cmp(&key_b.4))   // position 2
                .then_with(|| key_a.6.cmp(&key_b.6))            // sequence
        });
        if let Some((_key, max_vec)) = groups.first() {
            if max_vec.len() >= 2 {
                let consensus_record: &VariantRecord = max_vec
                    .iter()
                    .min_by_key(|vr| vr.get_read_id())
                    .unwrap();
                let read_ids: HashSet<ReadID> = max_vec.iter().map(|v| v.get_read_id()).collect();
                return (
                    consensus_record.get_graph_operation().clone(),
                    read_ids,
                    GraphOperationConsensusMethod::Modal
                );
            }
        }

        // Step 5. Otherwise, construct a consensus graph operation.
        // Positions 1 and 2 are median values.
        let chromosome_1: ReferenceChromosomeID = variant_records_.first().unwrap().get_chromosome_1();
        let chromosome_2: ReferenceChromosomeID = variant_records_.first().unwrap().get_chromosome_2();
        let strand_1: Strand = variant_records_.first().unwrap().get_strand_1().clone();
        let strand_2: Strand = variant_records_.first().unwrap().get_strand_2().clone();
        let operation_1: GraphOperationType = variant_records_.first().unwrap().get_operation_1().clone();
        let operation_2: GraphOperationType = variant_records_.first().unwrap().get_operation_2().clone();
        let variant_type: VariantType = variant_records_.first().unwrap().get_variant_type().clone();
        let mut position_1: Vec<ReferencePosition> = variant_records_.iter().map(|vr| vr.get_position_1()).collect();
        let mut position_2: Vec<ReferencePosition> = variant_records_.iter().map(|vr| vr.get_position_2()).collect();
        position_1.sort_unstable();
        position_2.sort_unstable();
        let position_1_median: ReferencePosition = position_1[(position_1.len() - 1) / 2];
        let position_2_median: ReferencePosition = position_2[(position_2.len() - 1) / 2];

        // Sequence is based on partial order alignment.
        let mut sequences: Vec<Box<str>> = variant_records_
            .iter()
            .map(|vr| vr.get_standardized_sequence())
            .filter(|sequence| !sequence.is_empty())
            .map(Box::from)
            .collect();
        let consensus_sequence: Box<str> = match sequences.len() {
            0 => "".into(),
            1 => sequences.pop().unwrap(),
            _ => perform_partial_order_alignment(
                &sequences,
                AlignmentMode::Global,
                match_score,
                mismatch_score,
                gap_open_score,
                gap_extend_score
            )
        };

        // The spellings aligned above are standardized, which is the forward orientation. A
        // graph operation holds its sequence in the orientation of strand 1, as a record does.
        let consensus_sequence: Box<str> = if strand_1 == Strand::Forward {
            consensus_sequence
        } else {
            reverse_complement(&consensus_sequence)
        };

        let consensus_graph_operation: GraphOperation = GraphOperation::new(
            chromosome_1,
            position_1_median,
            strand_1,
            operation_1,
            chromosome_2,
            position_2_median,
            strand_2,
            operation_2,
            consensus_sequence,
            variant_type
        );

        (consensus_graph_operation, HashSet::new(), GraphOperationConsensusMethod::Constructed)
    }

    pub fn get_alternate_allele_fraction(&self) -> f64 {
        assert!(self.total_depth > 0);
        let num_alt_allele_read_count: usize = self.get_read_ids().len();
        num_alt_allele_read_count as f64 / self.total_depth as f64
    }

    pub fn get_consensus_graph_operation(&self) -> &GraphOperation {
        &self.consensus_graph_operation
    }

    pub fn get_consensus_read_ids(&self) -> &HashSet<ReadID> {
        &self.consensus_read_ids
    }

    pub fn get_consensus_read_names(&self, read_names_map: &BiMap<ReadName,ReadID>) -> Vec<ReadName> {
        let mut read_ids: Vec<&ReadID> = self.get_consensus_read_ids().iter().collect();
        read_ids.sort_unstable();
        read_ids
            .into_iter()
            .map(|read_id| read_names_map.get_by_right(read_id).unwrap().clone())
            .collect()
    }

    pub fn get_consensus_method(&self) -> &GraphOperationConsensusMethod {
        &self.consensus_method
    }

    pub fn get_id(&self) -> VariantID {
        self.id
    }

    pub fn get_num_reads(&self) -> ReadSupport {
        self.get_read_ids().len() as ReadSupport
    }
    
    pub fn get_read_ids(&self) -> Vec<ReadID> {
        let mut read_ids: Vec<ReadID> = self.variant_records
            .iter()
            .map(|record| record.get_read_id())
            .collect();
        read_ids.sort_unstable();
        read_ids.dedup();
        read_ids
    }

    pub fn get_read_names(&self, read_names_map: &HashMap<ReadID,ReadName>) -> Vec<ReadName> {
        self.get_read_ids()
            .iter()
            .map(|read_id| read_names_map.get(read_id).unwrap().clone())
            .collect()
    }

    pub fn get_strand_1(&self) -> Strand {
        let has_forward: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.get_strand_1() == &Strand::Forward);
        let has_reverse: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.get_strand_1() == &Strand::Reverse);
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
            .any(|vr| vr.get_strand_2() == &Strand::Forward);
        let has_reverse: bool = self
            .variant_records
            .iter()
            .any(|vr| vr.get_strand_2() == &Strand::Reverse);
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

    pub fn get_variant_records(&self) -> &HashSet<VariantRecord> {
        &self.variant_records
    }

    pub fn into_variant_records(self) -> HashSet<VariantRecord> {
        self.variant_records
    }

    pub fn set_id(&mut self, id: VariantID) {
        self.id = id;
    }

    pub fn set_total_depth(&mut self, total_depth: i32) {
        self.total_depth = total_depth;
    }
}

impl Clone for VariantCall {
    fn clone(&self) -> Self {
        VariantCall {
            id: self.id,
            total_depth: self.total_depth,
            variant_records: self.variant_records.clone(),
            consensus_graph_operation: self.consensus_graph_operation.clone(),
            consensus_read_ids: self.consensus_read_ids.clone(),
            consensus_method: self.consensus_method.clone()
        }
    }
}


#[cfg(test)]
#[path = "../tests/variant/variant_call.rs"]
mod tests;