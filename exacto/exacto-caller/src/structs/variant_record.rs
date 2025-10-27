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
use serde::{Serialize, Deserialize};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Eq,Serialize,Deserialize)]
pub struct VariantRecord {
    pub read_id: usize,
    pub read_position_1: usize,
    pub read_position_2: usize,
    pub graph_operation: GraphOperation
}

impl PartialEq for VariantRecord {
    fn eq(&self, other: &Self) -> bool {
        if self.read_id == other.read_id &&
            self.read_position_1 == other.read_position_1 &&
            self.read_position_2 == other.read_position_2 &&
            self.graph_operation == other.graph_operation {
            true
        } else {
            false
        }
    }
}

impl Hash for VariantRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.read_id.hash(state);
        self.read_position_1.hash(state);
        self.read_position_2.hash(state);
        self.graph_operation.hash(state);
    }
}

impl VariantRecord {
    pub fn new(
        read_id: usize,
        read_start_position: usize,
        read_end_position: usize,
        sequence_operation: GraphOperation
    ) -> Self {
        Self {
            read_id,
            read_position_1: read_start_position,
            read_position_2: read_end_position,
            graph_operation: sequence_operation
        }
    }

    pub fn get_chromosome_1(&self) -> u16 {
        self.graph_operation.get_chromosome_1()
    }

    pub fn get_chromosome_2(&self) -> u16 {
        self.graph_operation.get_chromosome_2()
    }

    pub fn get_id(&self) -> Box<str> {
        format!("{}:{}-{}", self.read_id, self.read_position_1, self.read_position_2).into()
    }

    pub fn get_operation_1(&self) -> &GraphOperationType {
        self.graph_operation.get_operation_type_1()
    }

    pub fn get_operation_2(&self) -> &GraphOperationType {
        self.graph_operation.get_operation_type_2()
    }

    pub fn get_position_1(&self) -> usize {
        self.graph_operation.get_position_1()
    }

    pub fn get_position_2(&self) -> usize {
        self.graph_operation.get_position_2()
    }

    pub fn get_read_id(&self) -> usize {
        self.read_id
    }
    
    pub fn get_read_position_1(&self) -> usize {
        self.read_position_1
    }

    pub fn get_read_position_2(&self) -> usize {
        self.read_position_2
    }

    pub fn get_graph_operation(&self) -> &GraphOperation {
        &self.graph_operation
    }

    pub fn get_graph_operation_boxed_str(&self) -> Box<str> {
        self.graph_operation.as_boxed_str()
    }

    pub fn get_graph_operation_named_boxed_str(&self, chromosome_names_map: &BiMap<Box<str>,u16>) -> Box<str> {
        self.graph_operation.as_named_boxed_str(chromosome_names_map)
    }

    pub fn get_sequence(&self) -> &str {
        &*self.graph_operation.get_sequence()
    }

    pub fn get_sequence_length(&self) -> usize {
        self.graph_operation.get_sequence_length()
    }

    pub fn get_standardized_sequence(&self) -> String {
        self.graph_operation.get_standardized_sequence()
    }
    
    pub fn get_strand_1(&self) -> &Strand {
        self.graph_operation.get_strand_1()
    }

    pub fn get_strand_2(&self) -> &Strand {
        self.graph_operation.get_strand_2()
    }

    pub fn get_variant_size(&self) -> isize {
        self.graph_operation.get_variant_size()
    }

    pub fn get_variant_type(&self) -> &VariantType {
        self.graph_operation.get_variant_type()
    }
}

impl Clone for VariantRecord {
    fn clone(&self) -> Self {
        VariantRecord {
            read_id: self.read_id,
            read_position_1: self.read_position_1,
            read_position_2: self.read_position_2,
            graph_operation: self.graph_operation.clone()
        }
    }
}
