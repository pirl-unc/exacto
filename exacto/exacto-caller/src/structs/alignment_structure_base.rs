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


#[derive(Debug,Serialize,Deserialize)]
pub struct AlignmentStructureBase {
    read_position: u32,
    nucleotide: Nucleotide,
    base_quality: u8,
    kind: AlignmentStructureBaseKind,
    mapping_quality: Option<u16>,
    context: Option<AlignmentStructureBaseContext>,
    reference_chromosome_id: Option<u16>,
    reference_position: Option<u32>,
    reference_strand: Option<Strand>,
    reference_gene_id: Option<Box<str>>,
    reference_transcript_id: Option<Box<str>>,
    reference_exon_id: Option<Box<str>>
}

// API methods
impl AlignmentStructureBase {
    pub fn new(
        read_position: u32,
        nucleotide: Nucleotide,
        base_quality: u8
    ) -> Self {
        Self {
            read_position,
            nucleotide,
            base_quality,
            kind: AlignmentStructureBaseKind::Unaligned,
            // is_soft_clipped: false,
            // is_embedded_insertion: false,
            mapping_quality: None,
            context: None,
            reference_chromosome_id: None,
            reference_position: None,
            reference_strand: None,
            reference_gene_id: None,
            reference_transcript_id: None,
            reference_exon_id: None
        }
    }
    
    pub fn get_base_quality(&self) -> u8 {
        self.base_quality
    }

    pub fn get_context(&self) -> &Option<AlignmentStructureBaseContext> {
        &self.context
    }
    
    pub fn get_kind(&self) -> &AlignmentStructureBaseKind {
        &self.kind
    }

    pub fn get_mapping_quality(&self) -> &Option<u16> {
        &self.mapping_quality
    }

    pub fn get_nucleotide(&self) -> &Nucleotide {
        &self.nucleotide
    }

    pub fn get_read_position(&self) -> u32 {
        self.read_position
    }

    pub fn get_reference_exon_id(&self) -> &Option<Box<str>> {
        &self.reference_exon_id
    }

    pub fn get_reference_chromosome_id(&self) -> &Option<u16> {
        &self.reference_chromosome_id
    }

    pub fn get_reference_gene_id(&self) -> &Option<Box<str>> {
        &self.reference_gene_id
    }

    pub fn get_reference_position(&self) -> &Option<u32> {
        &self.reference_position
    }

    pub fn get_reference_strand(&self) -> &Option<Strand> {
        &self.reference_strand
    }

    pub fn get_reference_transcript_id(&self) -> &Option<Box<str>> {
        &self.reference_transcript_id
    }
    
    pub fn set_context(&mut self, context: AlignmentStructureBaseContext) {
        self.context = Some(context);
    }

    pub fn set_kind(&mut self, kind: AlignmentStructureBaseKind) {
        self.kind = kind;
    }
    
    pub fn set_mapping_quality(&mut self, value: u16) {
        self.mapping_quality = Some(value);
    }

    pub fn set_reference_exon_id(&mut self, value: &str) {
        self.reference_exon_id = Some(value.into());
    }

    pub fn set_reference_chromosome_id(&mut self, value: u16) {
        self.reference_chromosome_id = Some(value);
    }

    pub fn set_reference_gene_id(&mut self, value: &str) {
        self.reference_gene_id = Some(value.into());
    }
    
    pub fn set_reference_position(&mut self, value: u32) {
        self.reference_position = Some(value);
    }
    
    pub fn set_reference_strand(&mut self, value: Strand) {
        self.reference_strand = Some(value);
    }

    pub fn set_reference_transcript_id(&mut self, value: &str) {
        self.reference_transcript_id = Some(value.into());
    }
}

impl Clone for AlignmentStructureBase {
    fn clone(&self) -> Self {
        AlignmentStructureBase {
            read_position: self.read_position,
            nucleotide: self.nucleotide.clone(),
            base_quality: self.base_quality,
            kind: self.kind.clone(),
            mapping_quality: self.mapping_quality,
            context: self.context.clone(),
            reference_chromosome_id: self.reference_chromosome_id,
            reference_position: self.reference_position,
            reference_strand: self.reference_strand.clone(),
            reference_gene_id: self.reference_gene_id.clone(),
            reference_transcript_id: self.reference_transcript_id.clone(),
            reference_exon_id: self.reference_exon_id.clone()
        }
    }
}
