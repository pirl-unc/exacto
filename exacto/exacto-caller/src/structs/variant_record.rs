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
use serde::{Serialize, Deserialize};
use std::hash::{Hash, Hasher};

use crate::common::constants::*;
use crate::structs::sequence_operation::SequenceOperation;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct VariantRecord {
    pub read_id: usize,
    pub sequence_operation: SequenceOperation
}

impl Hash for VariantRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.read_id.hash(state);
        self.sequence_operation.hash(state);
    }
}

impl VariantRecord {
    pub fn new(
        read_id: usize,
        sequence_operation: SequenceOperation
    ) -> Self {
        Self {
            read_id,
            sequence_operation
        }
    }

    pub fn get_chromosome_1(&self) -> u16 {
        self.sequence_operation.chromosome_1
    }

    pub fn get_chromosome_2(&self) -> u16 {
        self.sequence_operation.chromosome_2
    }

    pub fn get_sequence_operation_boxed_str(&self) -> Box<str> {
        self.sequence_operation.as_boxed_str()
    }

    pub fn get_sequence_operation_named_boxed_str(&self, chromosome_names_map: &BiMap<Box<str>,u16>) -> Box<str> {
        self.sequence_operation.as_named_boxed_str(chromosome_names_map)
    }

    pub fn get_operation_1(&self) -> SequenceOperationTypes {
        self.sequence_operation.operation_1.clone()
    }

    pub fn get_operation_2(&self) -> SequenceOperationTypes {
        self.sequence_operation.operation_2.clone()
    }

    pub fn get_position_1(&self) -> u32 {
        self.sequence_operation.position_1
    }

    pub fn get_position_2(&self) -> u32 {
        self.sequence_operation.position_2
    }

    pub fn get_sequence(&self) -> &str {
        &*self.sequence_operation.sequence
    }

    pub fn get_sequence_length(&self) -> usize {
        self.sequence_operation.get_sequence_length()
    }

    pub fn get_strand_1(&self) -> Strands {
        self.sequence_operation.strand_1.clone()
    }

    pub fn get_strand_2(&self) -> Strands {
        self.sequence_operation.strand_2.clone()
    }

    pub fn get_variant_size(&self) -> isize {
        self.sequence_operation.get_variant_size()
    }

    pub fn get_variant_type(&self) -> SequenceOperationVariantTypes {
        self.sequence_operation.variant_type.clone()
    }
}

impl Clone for VariantRecord {
    fn clone(&self) -> Self {
        VariantRecord {
            read_id: self.read_id,
            sequence_operation: self.sequence_operation.clone()
        }
    }
}
