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
use serde::{Deserialize,Serialize};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct PrimaryStructureRecord {
    index: usize,
    record_type: PrimaryStructureRecordType,
    amino_acid: Option<Box<str>>,
    codon_index: Option<u8>,
    nucleotide: Option<Nucleotide>,
    transcript_model_id: usize,
    reference_transcript_ids: Vec<Box<str>>,
    transcript_structure_index: usize,
    read_start: usize,
    read_end: usize,
    net_variant_nucleotides_count: i32,
    frameshift_state: Option<FrameshiftState>,
    rna_variant_call_ids: Vec<usize>,
    dna_variant_call_ids: Vec<usize>
}

impl PrimaryStructureRecord {
    pub fn new(
        index: usize,
        record_type: PrimaryStructureRecordType,
        amino_acid: Option<Box<str>>,
        codon_index: Option<u8>,
        nucleotide: Option<Nucleotide>,
        transcript_model_id: usize,
        reference_transcript_ids: Vec<Box<str>>,
        transcript_structure_index: usize,
        read_start: usize,
        read_end: usize,
        net_variant_nucleotides_count: i32,
        frameshift_state: Option<FrameshiftState>,
        rna_variant_call_ids: Vec<usize>,
        dna_variant_call_ids: Vec<usize>
    ) -> Self {
        Self {
            index: index,
            record_type: record_type,
            amino_acid: amino_acid,
            codon_index: codon_index,
            nucleotide: nucleotide,
            transcript_model_id: transcript_model_id,
            reference_transcript_ids: reference_transcript_ids,
            transcript_structure_index: transcript_structure_index,
            read_start: read_start,
            read_end: read_end,
            net_variant_nucleotides_count: net_variant_nucleotides_count,
            frameshift_state: frameshift_state,
            rna_variant_call_ids: rna_variant_call_ids,
            dna_variant_call_ids: dna_variant_call_ids
        }
    }

    pub fn get_index(&self) -> usize {
        self.index
    }

    pub fn get_record_type(&self) -> &PrimaryStructureRecordType {
        &self.record_type
    }

    pub fn get_amino_acid(&self) -> &Option<Box<str>> {
        &self.amino_acid
    }

    pub fn get_codon_index(&self) -> &Option<u8> {
        &self.codon_index
    }

    pub fn get_nucleotide(&self) -> &Option<Nucleotide> {
        &self.nucleotide
    }

    pub fn get_transcript_model_id(&self) -> usize {
        self.transcript_model_id
    }

    pub fn get_reference_transcript_ids(&self) -> &Vec<Box<str>> {
        &self.reference_transcript_ids
    }

    pub fn get_transcript_structure_index(&self) -> usize {
        self.transcript_structure_index
    }

    pub fn get_read_start(&self) -> usize {
        self.read_start
    }

    pub fn get_read_end(&self) -> usize {
        self.read_end
    }

    pub fn get_net_variant_nucleotides_count(&self) -> i32 {
        self.net_variant_nucleotides_count
    }

    pub fn get_frameshift_state(&self) -> &Option<FrameshiftState> {
        &self.frameshift_state
    }

    pub fn get_rna_variant_call_ids(&self) -> &Vec<usize> {
        &self.rna_variant_call_ids
    }

    pub fn get_dna_variant_call_ids(&self) -> &Vec<usize> {
        &self.dna_variant_call_ids
    }

    pub fn set_amino_acid(&mut self, amino_acid: Option<Box<str>>) {
        self.amino_acid = amino_acid;
    }

    pub fn set_frameshift_state(&mut self, frameshift_state: Option<FrameshiftState>) {
        self.frameshift_state = frameshift_state;
    }
}

impl Clone for PrimaryStructureRecord {
    fn clone(&self) -> Self {
        PrimaryStructureRecord {
            index: self.index,
            record_type: self.record_type.clone(),
            amino_acid: self.amino_acid.clone(),
            codon_index: self.codon_index,
            nucleotide: self.nucleotide.clone(),
            transcript_model_id: self.transcript_model_id,
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            transcript_structure_index: self.transcript_structure_index,
            read_start: self.read_start,
            read_end: self.read_end,
            net_variant_nucleotides_count: self.net_variant_nucleotides_count,
            frameshift_state: self.frameshift_state.clone(),
            rna_variant_call_ids: self.rna_variant_call_ids.clone(),
            dna_variant_call_ids: self.dna_variant_call_ids.clone()
        }
    }
}
