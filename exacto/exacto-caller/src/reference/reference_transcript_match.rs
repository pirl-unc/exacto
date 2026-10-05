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


use exacto_core::prelude::{ReferenceGeneID, ReferenceGeneName, ReferenceTranscriptID};
use serde::{Deserialize, Serialize};
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct ReferenceTranscriptMatch {
    reference_gene_id: ReferenceGeneID,
    reference_gene_name: ReferenceGeneName,
    reference_transcript_id: ReferenceTranscriptID,
    splice_junction_matches: Vec<ReferenceTranscriptSpliceJunction>,
    num_overlapping_bases: u32,
    num_reference_only_bases: u32,
    num_query_only_bases: u32,
    num_terminal_query_only_bases: u32
}

impl PartialEq for ReferenceTranscriptMatch {
    fn eq(&self, other: &Self) -> bool {
        self.reference_gene_id == other.reference_gene_id 
            && self.reference_gene_name == other.reference_gene_name 
            && self.reference_transcript_id == other.reference_transcript_id 
            && self.splice_junction_matches == other.splice_junction_matches 
            && self.num_overlapping_bases == other.num_overlapping_bases 
            && self.num_reference_only_bases == other.num_reference_only_bases 
            && self.num_query_only_bases == other.num_query_only_bases 
            && self.num_terminal_query_only_bases == other.num_terminal_query_only_bases
    }
}

impl Eq for ReferenceTranscriptMatch {}

impl Hash for ReferenceTranscriptMatch {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_gene_id.hash(state);
        self.reference_gene_name.hash(state);
        self.reference_transcript_id.hash(state);
        self.splice_junction_matches.hash(state);
        self.num_overlapping_bases.hash(state);
        self.num_reference_only_bases.hash(state);
        self.num_query_only_bases.hash(state);
        self.num_terminal_query_only_bases.hash(state);
    }
}

impl ReferenceTranscriptMatch {
    pub fn new(
        reference_gene_id: &str,
        reference_gene_name: &str,
        reference_transcript_id: &str,
        splice_junction_matches: &Vec<ReferenceTranscriptSpliceJunction>,
        num_overlapping_bases: u32,
        num_reference_only_bases: u32,
        num_query_only_bases: u32,
        num_terminal_query_only_bases: u32
    ) -> Self {
        Self {
            reference_gene_id: reference_gene_id.into(),
            reference_gene_name: reference_gene_name.into(),
            reference_transcript_id: reference_transcript_id.into(),
            splice_junction_matches: splice_junction_matches.clone(),
            num_overlapping_bases: num_overlapping_bases,
            num_reference_only_bases: num_reference_only_bases,
            num_query_only_bases: num_query_only_bases,
            num_terminal_query_only_bases: num_terminal_query_only_bases
        }
    }
    
    pub fn get_reference_gene_id(&self) -> &str {
        &self.reference_gene_id
    }
    
    pub fn get_reference_gene_name(&self) -> &str {
        &self.reference_gene_name
    }
    
    pub fn get_reference_transcript_id(&self) -> &str {
        &self.reference_transcript_id
    }
    
    pub fn num_overlapping_bases(&self) -> u32 {
        self.num_overlapping_bases
    }
    
    pub fn num_reference_only_bases(&self) -> u32 {
        self.num_reference_only_bases
    }
    
    pub fn num_query_only_bases(&self) -> u32 {
        self.num_query_only_bases
    }

    pub fn num_terminal_query_only_bases(&self) -> u32 {
        self.num_terminal_query_only_bases
    }

    pub fn num_internal_query_only_bases(&self) -> u32 {
        self.num_query_only_bases.saturating_sub(self.num_terminal_query_only_bases)
    }

    pub fn num_splice_junction_matches(&self) -> usize {
        self.splice_junction_matches.len()
    }
}

impl Clone for ReferenceTranscriptMatch {
    fn clone(&self) -> Self {
        ReferenceTranscriptMatch {
            reference_gene_id: self.reference_gene_id.clone(),
            reference_gene_name: self.reference_gene_name.clone(),
            reference_transcript_id: self.reference_transcript_id.clone(),
            splice_junction_matches: self.splice_junction_matches.clone(),
            num_overlapping_bases: self.num_overlapping_bases,
            num_reference_only_bases: self.num_reference_only_bases,
            num_query_only_bases: self.num_query_only_bases,
            num_terminal_query_only_bases: self.num_terminal_query_only_bases
        }
    }
}
