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


use serde::{Deserialize, Serialize};


#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct IntegratedVariantRecord {
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub rna_variant_id: u32,
    pub dna_variant_id: u32,
    pub distance: u32,
    pub rna_variant_position: Box<str>,
    pub dna_variant_position: Box<str>
}
