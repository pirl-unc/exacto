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


use exacto_core::prelude::*;
use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct ReferenceBase {
    pub reference_chromosome_id: u16,
    pub reference_position: usize,
    
    /// Reference nucleotide on the reference (forward or reverse) strand.
    pub reference_nucleotide: Nucleotide,

    pub reference_strand: Strand,
    
    pub reference_gene_id: Option<Box<str>>,
    pub reference_transcript_id: Option<Box<str>>,
    pub reference_exon_id: Option<Box<str>>,
}

/// API methods
impl ReferenceBase {
    pub fn new(
        reference_chromosome_id: u16,
        reference_position: usize,
        reference_nucleotide: Nucleotide,
        reference_strand: Strand,
        reference_gene_id: Option<Box<str>>,
        reference_transcript_id: Option<Box<str>>,
        reference_exon_id: Option<Box<str>>
    ) -> Self {
        Self {
            reference_chromosome_id: reference_chromosome_id,
            reference_position: reference_position,
            reference_nucleotide: reference_nucleotide,
            reference_strand: reference_strand,
            reference_gene_id: reference_gene_id,
            reference_transcript_id: reference_transcript_id,
            reference_exon_id: reference_exon_id
        }
    }
}

impl Clone for ReferenceBase {
    fn clone(&self) -> Self {
        ReferenceBase {
            reference_chromosome_id: self.reference_chromosome_id,
            reference_position: self.reference_position,
            reference_nucleotide: self.reference_nucleotide.clone(),
            reference_strand: self.reference_strand.clone(),
            reference_gene_id: self.reference_gene_id.clone(),
            reference_transcript_id: self.reference_transcript_id.clone(),
            reference_exon_id: self.reference_exon_id.clone()
        }
    }
}
