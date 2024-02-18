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
pub struct VariantCall {
    pub id: String,
    pub sample_id: String,
    pub nucleic_acid: String,
    pub chromosome_1: String,
    pub position_1: u32,
    pub chromosome_2: String,
    pub position_2: u32,
    pub variant_type: String,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub variant_size: u32,
    pub alternate_allele_read_ids: Vec<String>
}

impl VariantCall {
    pub fn new(
        id: String,
        sample_id: String,
        nucleic_acid: String,
        chromosome_1: String,
        position_1: u32,
        chromosome_2: String,
        position_2: u32,
        variant_type: String,
        reference_allele: String,
        alternate_allele: String,
        variant_size: u32
    ) -> Self {
        Self {
            id: id,
            sample_id: sample_id,
            nucleic_acid: nucleic_acid,
            chromosome_1: chromosome_1,
            position_1: position_1,
            chromosome_2: chromosome_2,
            position_2: position_2,
            variant_type: variant_type,
            reference_allele: reference_allele,
            alternate_allele: alternate_allele,
            variant_size: variant_size,
            alternate_allele_read_ids: Vec::new()
        }
    }

    pub fn add_alternate_allele_read_id(&mut self, alternate_allele_read_id: String) {
        self.alternate_allele_read_ids.push(alternate_allele_read_id);
    }
}

impl Clone for VariantCall {
    fn clone(&self) -> Self {
        VariantCall {
            id: self.id.to_string(),
            sample_id: self.sample_id.to_string(),
            nucleic_acid: self.nucleic_acid.to_string(),
            chromosome_1: self.chromosome_1.to_string(),
            position_1: self.position_1,
            chromosome_2: self.chromosome_2.to_string(),
            position_2: self.position_2,
            variant_type: self.variant_type.to_string(),
            reference_allele: self.reference_allele.to_string(),
            alternate_allele: self.alternate_allele.to_string(),
            variant_size: self.variant_size,
            alternate_allele_read_ids: self.alternate_allele_read_ids.clone()
        }
    }
}
