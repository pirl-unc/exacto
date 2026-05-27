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

use crate::prelude::*;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct IntegratedVariant {
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub reference_gene_names: Vec<Box<str>>,
    pub reference_transcript_ids: Vec<Box<str>>,
    pub rna_variant_id: u32,
    
    /// A map between DNA variant call ID and IntegratedVariantDistance object
    pub dna_variant_ids: HashMap<u32, IntegratedVariantDistance>
}

impl Hash for IntegratedVariant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.assembled_transcript_name.hash(state);
        self.transcript_model_id.hash(state);
        self.rna_variant_id.hash(state);
        self.reference_gene_names.hash(state);
        self.reference_transcript_ids.hash(state);
        for key in self.dna_variant_ids.keys() {
            key.hash(state);
        }
    }
}

impl IntegratedVariant {
    pub fn new(
        assembled_transcript_name: Box<str>,
        transcript_model_id: u32,
        reference_gene_names: &Vec<Box<str>>,
        reference_transcript_ids: &Vec<Box<str>>,
        rna_variant_id: u32
    ) -> Self {
        Self {
            assembled_transcript_name,
            transcript_model_id,
            reference_gene_names: reference_gene_names.clone(),
            reference_transcript_ids: reference_transcript_ids.clone(),
            rna_variant_id,
            dna_variant_ids: HashMap::new()
        }
    }
    
    pub fn add_dna_variant_id(
        &mut self, 
        dna_variant_id: u32, 
        distance: IntegratedVariantDistance
    ) {
        assert_eq!(self.dna_variant_ids.contains_key(&dna_variant_id), false);
        self.dna_variant_ids.insert(dna_variant_id, distance);
    }
}

impl Clone for IntegratedVariant {
    fn clone(&self) -> Self {
        IntegratedVariant {
            assembled_transcript_name: self.assembled_transcript_name.clone(),
            transcript_model_id: self.transcript_model_id,
            reference_gene_names: self.reference_gene_names.clone(),
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            rna_variant_id: self.rna_variant_id,
            dna_variant_ids: self.dna_variant_ids.clone()
        }
    }
}
