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


use exacto_core::prelude::ReadID;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::prelude::*;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RNAReadVariantRecordSet {
    read_id: ReadID,
    variant_records: HashSet<VariantRecord>
}

impl RNAReadVariantRecordSet {
    pub fn new(read_id: ReadID) -> Self {
        RNAReadVariantRecordSet {
            read_id: read_id,
            variant_records: HashSet::new()
        }
    }

    pub fn add_variant_record(&mut self, record: VariantRecord) {
        self.variant_records.insert(record);
    }

    pub fn get_read_id(&self) -> ReadID {
        self.read_id
    }

    pub fn get_variant_records(&self) -> &HashSet<VariantRecord> {
        &self.variant_records
    }
}