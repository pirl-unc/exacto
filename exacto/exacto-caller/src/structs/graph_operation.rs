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


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct GraphOperation {
    chromosome_1: u16,
    position_1: usize,
    strand_1: Strand,
    operation_type_1: GraphOperationType,
    chromosome_2: u16,
    position_2: usize,
    strand_2: Strand,
    operation_type_2: GraphOperationType,
    sequence: Box<str>,
    variant_type: VariantType
}

impl Hash for GraphOperation {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome_1.hash(state);
        self.position_1.hash(state);
        self.strand_1.hash(state);
        self.operation_type_1.hash(state);
        self.chromosome_2.hash(state);
        self.position_2.hash(state);
        self.strand_2.hash(state);
        self.operation_type_1.hash(state);
        self.sequence.hash(state);
        self.variant_type.hash(state);
    }
}

/// Constructor and GET functions
impl GraphOperation {
    pub fn new(
        chromosome_1: u16,
        position_1: usize,
        strand_1: Strand,
        operation_type_1: GraphOperationType,
        chromosome_2: u16,
        position_2: usize,
        strand_2: Strand,
        operation_type_2: GraphOperationType,
        sequence: Box<str>,
        variant_type: VariantType
    ) -> Self {
        let unstandardized: GraphOperation = Self {
            chromosome_1,
            position_1,
            strand_1,
            operation_type_1,
            chromosome_2,
            position_2,
            strand_2,
            operation_type_2,
            sequence,
            variant_type,
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
        s.push_str(":");
        s.push_str(self.variant_type.as_str());
        s.into_boxed_str()
    }

    pub fn as_named_boxed_str(&self, chromosome_names_map: &BiMap<Box<str>,u16>) -> Box<str> {
        let mut s: String = String::new();
        s.push_str(chromosome_names_map.get_by_right(&self.chromosome_1).unwrap());
        s.push_str(":");
        s.push_str(self.position_1.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_1.as_str());
        s.push_str(":");
        s.push_str(self.operation_type_1.as_str());
        s.push_str(":");
        s.push_str(chromosome_names_map.get_by_right(&self.chromosome_2).unwrap());
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
        s.push_str(":");
        s.push_str(self.variant_type.as_str());
        s.into_boxed_str()
    }

    pub fn get_chromosome_1(&self) -> u16 {
        self.chromosome_1
    }
    
    pub fn get_chromosome_2(&self) -> u16 {
        self.chromosome_2
    }
    
    pub fn get_position_1(&self) -> usize {
        self.position_1
    }
    
    pub fn get_position_2(&self) -> usize {
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

    pub fn get_standardized_sequence(&self) -> String {
        if self.strand_1 == Strand::Forward {
            self.sequence.to_string().to_uppercase()
        } else {
            reverse_complement(&self.sequence).to_string().to_uppercase()
        }
    }

    pub fn get_variant_size(&self) -> isize {
        match self.variant_type {
            VariantType::Breakpoint => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            VariantType::CircularRNA => {
                -1
            },
            VariantType::CrypticExon => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            VariantType::Deletion => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            VariantType::ExonTruncation => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            }
            VariantType::FusionGene => {
                -1
            },
            VariantType::Insertion => {
                self.get_sequence_length() as isize
            },
            VariantType::IntronRetention => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            VariantType::MultiNucleotideVariant => {
                self.get_sequence_length() as isize
            },
            VariantType::SingleNucleotideVariant => {
                1
            },
            VariantType::Translocation => {
                -1
            },
            VariantType::UTRExtension => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            }
        }
    }
    
    pub fn get_variant_type(&self) -> &VariantType {
        &self.variant_type
    }
}

/// Helper methods
impl GraphOperation {
    fn standardize(op: GraphOperation) -> GraphOperation {
        let should_swap: bool = if op.chromosome_1 != op.chromosome_2 {
            op.chromosome_1 > op.chromosome_2
        } else {
            op.position_1 > op.position_2
        };

        if should_swap {
            GraphOperation {
                chromosome_1: op.chromosome_2,
                position_1: op.position_2,
                strand_1: op.strand_2,
                operation_type_1: op.operation_type_2,
                chromosome_2: op.chromosome_1,
                position_2: op.position_1,
                strand_2: op.strand_1,
                operation_type_2: op.operation_type_1,
                sequence: op.sequence,
                variant_type: op.variant_type
            }
        } else {
            op
        }
    }
}

impl Clone for GraphOperation {
    fn clone(&self) -> Self {
        GraphOperation {
            chromosome_1: self.chromosome_1,
            position_1: self.position_1,
            strand_1: self.strand_1.clone(),
            operation_type_1: self.operation_type_1.clone(),
            chromosome_2: self.chromosome_2,
            position_2: self.position_2,
            strand_2: self.strand_2.clone(),
            operation_type_2: self.operation_type_2.clone(),
            sequence: self.sequence.clone(),
            variant_type: self.variant_type.clone()
        }
    }
}
