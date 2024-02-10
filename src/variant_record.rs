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


extern crate serde;
use serde::{Serialize, Deserialize};


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantRecord {
    pub id: String,
    pub read_id: String,
    pub chromosome_1: String,
    pub position_1: isize,
    pub chromosome_2: String,
    pub position_2: isize,
    pub variant_type: String,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub variant_size: isize
}

impl VariantRecord {
    pub fn new(
        id: String,
        read_id: String,
        chromosome_1: String,
        position_1: isize,
        chromosome_2: String,
        position_2: isize,
        variant_type: String,
        reference_allele: String,
        alternate_allele: String,
        variant_size: isize
    ) -> Self {
        Self {
            id: id,
            read_id: read_id,
            chromosome_1: chromosome_1,
            position_1: position_1,
            chromosome_2: chromosome_2,
            position_2: position_2,
            variant_type: variant_type,
            reference_allele: reference_allele,
            alternate_allele: alternate_allele,
            variant_size: variant_size
        }
    }
}

impl Clone for VariantRecord {
    fn clone(&self) -> Self {
        VariantRecord {
            id: self.id.clone(),
            read_id: self.read_id.clone(),
            chromosome_1: self.chromosome_1.clone(),
            position_1: self.position_1,
            chromosome_2: self.chromosome_2.clone(),
            position_2: self.position_2,
            variant_type: self.variant_type.clone(),
            reference_allele: self.reference_allele.clone(),
            alternate_allele: self.alternate_allele.clone(),
            variant_size: self.variant_size
        }
    }
}
