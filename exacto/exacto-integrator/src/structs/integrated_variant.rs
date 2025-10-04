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
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::structs::integrated_variant_distance::IntegratedVariantDistance;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct IntegratedVariant {
    pub transcript_model_id: usize,
    pub reference_transcript_ids: Vec<Box<str>>,
    pub rna_variant_call_id: usize,
    
    /// A map between DNA variant call ID and IntegratedVariantDistance object
    pub dna_variant_call_ids: HashMap<usize, IntegratedVariantDistance>
}

impl Hash for IntegratedVariant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.transcript_model_id.hash(state);
        self.rna_variant_call_id.hash(state);
        self.reference_transcript_ids.hash(state);
        for key in self.dna_variant_call_ids.keys() {
            key.hash(state);
        }
    }
}

impl IntegratedVariant {
    pub fn new(transcript_model_id: usize, reference_transcript_ids: &Vec<Box<str>>, rna_variant_call_id: usize) -> Self {
        Self {
            transcript_model_id: transcript_model_id,
            reference_transcript_ids: reference_transcript_ids.clone(),
            rna_variant_call_id: rna_variant_call_id,
            dna_variant_call_ids: HashMap::new()
        }
    }
    
    pub fn add_dna_variant_call_id(&mut self, dna_variant_call_id: usize, distance: IntegratedVariantDistance) {
        assert_eq!(self.dna_variant_call_ids.contains_key(&dna_variant_call_id), false);
        self.dna_variant_call_ids.insert(dna_variant_call_id, distance);
    }
}

impl Clone for IntegratedVariant {
    fn clone(&self) -> Self {
        IntegratedVariant {
            transcript_model_id: self.transcript_model_id,
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            rna_variant_call_id: self.rna_variant_call_id,
            dna_variant_call_ids: self.dna_variant_call_ids.clone()
        }
    }
}
