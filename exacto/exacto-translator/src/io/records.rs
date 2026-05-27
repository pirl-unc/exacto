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
use std::collections::HashSet;


#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptSupportRecord {
    pub assembled_transcript_name: Box<str>,
    pub sequence: Box<str>,
    pub read_names: Box<str>
}


#[derive(Debug, Serialize)]
pub struct NucleotideRecord {
    pub primary_structure_id: u32,
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub amino_acid_index: u32,
    pub amino_acid: Box<str>,
    pub codon_index: u8,
    pub nucleotide: Box<str>,
    pub is_amino_acid_variant: bool,
    pub is_nucleotide_variant: bool,
    pub transcript_read_position: u32,
    pub transcript_structure_index: Option<u32>,
    pub rna_variant_id: Option<u32>,
    pub rna_variant: Option<Box<str>>,
    pub dna_variant_ids: Option<HashSet<u32>>,
    pub dna_variant: Option<Box<str>>,
    pub preceding_event_rna_variant_id: Option<u32>,
    pub preceding_event_rna_variant: Option<Box<str>>,
    pub preceding_event_dna_variant_ids: Option<HashSet<u32>>,
    pub preceding_event_dna_variant: Option<Box<str>>
}


#[derive(Debug, Serialize)]
pub struct PrimaryStructureRecord {
    pub primary_structure_id: usize,
    pub amino_acid_sequence: String,
    pub amino_acid_sequence_length: usize,
    pub num_mutant_amino_acids: u32,
    pub assembled_transcript_name: Box<str>,
    pub assembled_transcript_sequence: Box<str>,
    pub assembled_transcript_sequence_length: usize,
    pub orf_start: u32,
    pub orf_end: u32,
    pub transcript_model_id: u32,
    pub reference_gene_names: String,
    pub reference_transcript_ids: String,
    pub mutant_amino_acid_intervals: String,
    pub rna_variant_ids: String,
    pub rna_variants: String,
    pub dna_variant_ids: String,
    pub dna_variants: String,
    pub rna_read_names: Box<str>,
    pub num_rna_read_names: usize,
    pub dna_variant_read_names: Box<str>,
    pub transcript_read_position: Box<str>
}