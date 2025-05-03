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


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct SequenceOperation {
    pub chromosome_1: u16,
    pub position_1: u32,
    pub strand_1: Strand,
    pub operation_1: SequenceOperationType,
    pub chromosome_2: u16,
    pub position_2: u32,
    pub strand_2: Strand,
    pub operation_2: SequenceOperationType,
    pub sequence: Box<str>,
    pub variant_type: SequenceOperationVariantType
}

impl Hash for SequenceOperation {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome_1.hash(state);
        self.position_1.hash(state);
        self.strand_1.hash(state);
        self.operation_1.hash(state);
        self.chromosome_2.hash(state);
        self.position_2.hash(state);
        self.strand_2.hash(state);
        self.operation_1.hash(state);
        self.sequence.hash(state);
        self.variant_type.hash(state);
    }
}

impl SequenceOperation {
    pub fn new(
        chromosome_1: u16,
        position_1: u32,
        strand_1: Strand,
        operation_1: SequenceOperationType,
        chromosome_2: u16,
        position_2: u32,
        strand_2: Strand,
        operation_2: SequenceOperationType,
        sequence: Box<str>,
        variant_type: SequenceOperationVariantType
    ) -> Self {
        Self {
            chromosome_1,
            position_1,
            strand_1,
            operation_1,
            chromosome_2,
            position_2,
            strand_2,
            operation_2,
            sequence,
            variant_type
        }
    }

    pub fn as_boxed_str(&self) -> Box<str> {
        let mut s: String = String::new();
        s.push_str(self.chromosome_1.to_string().as_str());
        s.push_str(":");
        s.push_str(self.position_1.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_1.as_str());
        s.push_str(":");
        s.push_str(self.operation_1.as_str());
        s.push_str(":");
        s.push_str(self.chromosome_2.to_string().as_str());
        s.push_str(":");
        s.push_str(self.position_2.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_2.as_str());
        s.push_str(":");
        s.push_str(self.operation_2.as_str());
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
        s.push_str(self.operation_1.as_str());
        s.push_str(":");
        s.push_str(chromosome_names_map.get_by_right(&self.chromosome_2).unwrap());
        s.push_str(":");
        s.push_str(self.position_2.to_string().as_str());
        s.push_str(":");
        s.push_str(self.strand_2.as_str());
        s.push_str(":");
        s.push_str(self.operation_2.as_str());
        s.push_str(":");
        s.push_str(&*self.sequence);
        s.push_str(":");
        s.push_str(self.get_sequence_length().to_string().as_str());
        s.push_str(":");
        s.push_str(self.variant_type.as_str());
        s.into_boxed_str()
    }

    pub fn get_sequence_length(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_variant_size(&self) -> isize {
        match self.variant_type {
            SequenceOperationVariantType::Alternative3PrimeSpliceSite => {
                -1
            },
            SequenceOperationVariantType::Alternative5PrimeSpliceSite => {
                -1
            },
            SequenceOperationVariantType::Breakpoint => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            SequenceOperationVariantType::CrypticExon => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            SequenceOperationVariantType::Deletion => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            SequenceOperationVariantType::ExonSkipping => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            SequenceOperationVariantType::FusionGene => {
                -1
            },
            SequenceOperationVariantType::Insertion => {
                self.get_sequence_length() as isize
            },
            SequenceOperationVariantType::IntronRetention => {
                (self.position_2.abs_diff(self.position_1) as isize) - 1
            },
            SequenceOperationVariantType::MultiNucleotideVariant => {
                self.get_sequence_length() as isize
            },
            SequenceOperationVariantType::SingleNucleotideVariant => {
                1
            },
            SequenceOperationVariantType::Translocation => {
                -1
            }
        }
    }
}

impl Clone for SequenceOperation {
    fn clone(&self) -> Self {
        SequenceOperation {
            chromosome_1: self.chromosome_1,
            position_1: self.position_1,
            strand_1: self.strand_1.clone(),
            operation_1: self.operation_1.clone(),
            chromosome_2: self.chromosome_2,
            position_2: self.position_2,
            strand_2: self.strand_2.clone(),
            operation_2: self.operation_2.clone(),
            sequence: self.sequence.clone(),
            variant_type: self.variant_type.clone()
        }
    }
}
