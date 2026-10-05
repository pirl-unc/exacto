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
use std::collections::BTreeMap;

use crate::prelude::*;


#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub struct RNAVariantIntegration {
    pub assembled_transcript_name: Box<str>,
    pub reference_gene_names: Vec<Box<str>>,
    pub reference_transcript_ids: Vec<Box<str>>,
    pub rna_variant_id: u32,
    
    /// A map between DNA variant call ID and DNAVariantMatch object, in ascending ID order.
    pub dna_variant_ids: BTreeMap<u32, DNAVariantMatch>
}

impl RNAVariantIntegration {
    pub fn new(
        assembled_transcript_name: Box<str>,
        reference_gene_names: Vec<Box<str>>,
        reference_transcript_ids: Vec<Box<str>>,
        rna_variant_id: u32
    ) -> Self {
        Self {
            assembled_transcript_name: assembled_transcript_name,
            reference_gene_names: reference_gene_names,
            reference_transcript_ids: reference_transcript_ids,
            rna_variant_id: rna_variant_id,
            dna_variant_ids: BTreeMap::new()
        }
    }
    
    pub fn add_dna_variant_id(
        &mut self, 
        dna_variant_id: u32,
        dna_variant_match: DNAVariantMatch
    ) {
        assert_eq!(self.dna_variant_ids.contains_key(&dna_variant_id), false);
        self.dna_variant_ids.insert(dna_variant_id, dna_variant_match);
    }
}


#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct DNAVariantMatch {
    distance: u32,
    rna_variant_position_index: VariantPosition,
    dna_variant_position_index: VariantPosition
}

impl DNAVariantMatch {
    pub fn new(
        distance: u32,
        rna_variant_position_used: VariantPosition,
        dna_variant_position_used: VariantPosition
    ) -> Self {
        Self {
            distance: distance,
            rna_variant_position_index: rna_variant_position_used,
            dna_variant_position_index: dna_variant_position_used
        }
    }
    
    pub fn get_distance(&self) -> u32 {
        self.distance
    }
    
    pub fn get_rna_variant_position_index(&self) -> &VariantPosition {
        &self.rna_variant_position_index
    }
    
    pub fn get_dna_variant_position_index(&self) -> &VariantPosition {
        &self.dna_variant_position_index
    }
}
