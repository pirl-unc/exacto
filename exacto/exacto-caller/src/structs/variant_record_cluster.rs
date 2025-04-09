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


use std::sync::Arc;

use crate::prelude::*;


#[derive(Debug,Clone)]
pub struct VariantRecordCluster {
    pub variant_records: Vec<Arc<VariantRecord>>,
    pub chromosome_1: u16,
    pub chromosome_2: u16,
    pub min_position_1: u32,
    pub max_position_1: u32,
    pub min_position_2: u32,
    pub max_position_2: u32
}

impl VariantRecordCluster {
    pub fn new(
        chromosome_1: u16,
        chromosome_2: u16,
        min_position_1: u32,
        max_position_1: u32,
        min_position_2: u32,
        max_position_2: u32
    ) -> Self {
        Self {
            variant_records: Vec::new(),
            chromosome_1: chromosome_1,
            chromosome_2: chromosome_2,
            min_position_1: min_position_1,
            max_position_1: max_position_1,
            min_position_2: min_position_2,
            max_position_2: max_position_2
        }
    }

    pub fn add_variant_record(&mut self, variant_record: Arc<VariantRecord>) {
        if self.variant_records.is_empty() {
            self.variant_records.push(variant_record);
        } else {
            if variant_record.get_variant_type() == self.variant_records[0].get_variant_type() {
                self.variant_records.push(variant_record);
            } else {
                panic!("The variant record's type must be the same as the existing variant records.");
            }
        }
    }

    pub fn get_variant_type(&self) -> SequenceOperationVariantTypes {
        self.variant_records[0].get_variant_type()
    }

    pub fn overlaps_position_1(&self, other: &VariantRecordCluster) -> bool {
        if self.chromosome_1 == other.chromosome_1 &&
            other.max_position_1 >= self.min_position_1 &&
            other.min_position_1 <= self.max_position_1 {
            true
        } else {
            false
        }
    }

    pub fn overlaps_position_2(&self, other: &VariantRecordCluster) -> bool {
        if self.chromosome_2 == other.chromosome_2 &&
            other.max_position_2 >= self.min_position_2 &&
            other.min_position_2 <= self.max_position_2 {
            true
        } else {
            false
        }
    }
}
