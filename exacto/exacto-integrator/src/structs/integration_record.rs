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


use serde::{Serialize, Deserialize};


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct IntegrationRecord {
    pub rna_variant_id: usize,
    pub dna_variant_id: usize

}

impl Hash for IntegrationRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.rna_variant_id.hash(state);
        self.dna_variant_id.hash(state);
    }
}

impl IntegrationRecord {
    pub fn new(
        rna_variant_id: usize,
        dna_variant_id: usize
    ) -> Self {
        Self {
            rna_variant_id,
            dna_variant_id
        }
    }
}

impl Clone for IntegrationRecord {
    fn clone(&self) -> Self {
        IntegrationRecord {
            rna_variant_id: self.read_id,
            dna_variant_id: self.dna_variant_id
        }
    }
}
