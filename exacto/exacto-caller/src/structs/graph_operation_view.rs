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
use serde::{Serialize, Deserialize};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::prelude::*;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct GraphOperationView {
    chromosome_1: Arc<str>,
    position_1: u32,
    strand_1: Strand,
    operation_type_1: GraphOperationType,
    chromosome_2: Arc<str>,
    position_2: u32,
    strand_2: Strand,
    operation_type_2: GraphOperationType,
    sequence: Arc<str>,
    num_cycles: Option<u16>
}

impl Hash for GraphOperationView {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome_1.hash(state);
        self.position_1.hash(state);
        self.strand_1.hash(state);
        self.operation_type_1.hash(state);
        self.chromosome_2.hash(state);
        self.position_2.hash(state);
        self.strand_2.hash(state);
        self.operation_type_2.hash(state);
        self.sequence.hash(state);
        self.num_cycles.hash(state);
    }
}

impl GraphOperationView {
    pub fn new(
        chromosome_1: &str,
        position_1: u32,
        operation_type_1: &GraphOperationType,
        strand_1: &Strand,
        chromosome_2: &str,
        position_2: u32,
        operation_type_2: &GraphOperationType,
        strand_2: &Strand,
        sequence: &str,
        num_cycles: Option<u16>
    ) -> Self {
        if *operation_type_1 == GraphOperationType::Upstream &&
            *operation_type_2 == GraphOperationType::Downstream &&
            chromosome_1 == chromosome_2 &&
            position_1 < position_2 {
            assert!(
                num_cycles.is_some() == true,
                "A cycle creating operation must have a number of cycles: {}:{}:{}:{}:{}:{}:{}:{}:{}:{:?}",
                chromosome_1,
                position_1,
                operation_type_1.as_str(),
                strand_1.as_str(),
                chromosome_2,
                position_2,
                operation_type_2.as_str(),
                strand_2.as_str(),
                sequence,
                num_cycles
            );
        }

        let unstandardized: GraphOperationView = Self {
            chromosome_1: chromosome_1.into(),
            position_1: position_1,
            operation_type_1: operation_type_1.clone(),
            strand_1: strand_1.clone(),
            chromosome_2: chromosome_2.into(),
            position_2: position_2,
            operation_type_2: operation_type_2.clone(),
            strand_2: strand_2.clone(),
            sequence: sequence.into(),
            num_cycles: num_cycles
        };
        Self::standardize(unstandardized)
    }

    pub fn as_boxed_str(&self) -> Box<str> {
        let mut s: String = String::new();
        s.push_str(self.chromosome_1.to_string().as_str());
        s.push_str(":");
        s.push_str(self.position_1.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_1.as_str());
        s.push_str(":");
        s.push_str(self.operation_type_1.as_str());
        s.push_str(":");
        s.push_str(self.chromosome_2.to_string().as_str());
        s.push_str(":");
        s.push_str(self.position_2.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_2.as_str());
        s.push_str(":");
        s.push_str(self.operation_type_2.as_str());
        s.push_str(":");
        s.push_str(&*self.sequence);
        s.push_str(":");
        s.push_str(self.get_sequence_length().to_string().as_str());
        s.into_boxed_str()
    }

    pub fn get_chromosome_1(&self) -> &str {
        &*self.chromosome_1
    }

    pub fn get_chromosome_2(&self) -> &str {
        &*self.chromosome_2
    }

    pub fn get_position_1(&self) -> u32 {
        self.position_1
    }

    pub fn get_position_2(&self) -> u32 {
        self.position_2
    }

    pub fn get_strand_1(&self) -> &Strand {
        &self.strand_1
    }

    pub fn get_strand_2(&self) -> &Strand {
        &self.strand_2
    }

    pub fn get_operation_type_1(&self) -> &GraphOperationType {
        &self.operation_type_1
    }

    pub fn get_operation_type_2(&self) -> &GraphOperationType {
        &self.operation_type_2
    }

    pub fn get_sequence(&self) -> &str {
        &*self.sequence
    }

    pub fn get_sequence_length(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_standardized_sequence(&self) -> Box<str> {
        if self.strand_1 == Strand::Forward {
            self.sequence.to_string().to_uppercase().into()
        } else {
            reverse_complement(&self.sequence).to_string().to_uppercase().into()
        }
    }

    pub fn get_num_cycles(&self) -> &Option<u16> {
        &self.num_cycles
    }

    pub fn get_variant_type(&self) -> VariantType {
        if self.chromosome_1 == self.chromosome_2 {
            if self.operation_type_1 == GraphOperationType::Downstream &&
                self.operation_type_2 == GraphOperationType::Upstream {
                if self.sequence.len() == 0 {
                    VariantType::Deletion
                } else if self.sequence.len() == 1 {
                    if self.position_1.abs_diff(self.position_2) == 2 {
                        VariantType::SingleNucleotideVariant
                    } else {
                        VariantType::Insertion
                    }
                } else {
                    if self.position_1.abs_diff(self.position_2) == self.sequence.len() as u32 + 1 {
                        VariantType::MultiNucleotideVariant
                    } else {
                        VariantType::Insertion
                    }
                }
            } else if self.operation_type_1 == GraphOperationType::Upstream && self.operation_type_2 == GraphOperationType::Downstream {
                VariantType::Breakpoint
            } else if self.operation_type_1 == GraphOperationType::Downstream && self.operation_type_2 == GraphOperationType::Downstream {
                VariantType::Breakpoint
            } else if self.operation_type_1 == GraphOperationType::Upstream && self.operation_type_2 == GraphOperationType::Upstream {
                VariantType::Breakpoint
            } else {
                panic!("Unsupported operation types: {:?} {:?}", self.operation_type_1, self.operation_type_2)
            }
        } else {
            VariantType::Translocation
        }
    }
}

/// Helper methods
impl GraphOperationView {
    fn standardize(op: GraphOperationView) -> GraphOperationView {
        let should_swap: bool = if op.chromosome_1 != op.chromosome_2 {
            op.chromosome_1 > op.chromosome_2
        } else {
            op.position_1 > op.position_2
        };

        if should_swap {
            GraphOperationView {
                chromosome_1: op.chromosome_2,
                position_1: op.position_2,
                operation_type_1: op.operation_type_2,
                strand_1: op.strand_2,
                chromosome_2: op.chromosome_1,
                position_2: op.position_1,
                operation_type_2: op.operation_type_1,
                strand_2: op.strand_1,
                sequence: op.sequence,
                num_cycles: op.num_cycles
            }
        } else {
            op
        }
    }
}

impl Clone for GraphOperationView {
    fn clone(&self) -> Self {
        GraphOperationView {
            chromosome_1: self.chromosome_1.clone(),
            position_1: self.position_1,
            operation_type_1: self.operation_type_1.clone(),
            strand_1: self.strand_1.clone(),
            chromosome_2: self.chromosome_2.clone(),
            position_2: self.position_2,
            operation_type_2: self.operation_type_2.clone(),
            strand_2: self.strand_2.clone(),
            sequence: self.sequence.clone(),
            num_cycles: self.num_cycles
        }
    }
}
