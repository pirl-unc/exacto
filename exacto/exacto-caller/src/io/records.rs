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


#[derive(Debug, Serialize, Clone)]
pub struct AssembledTranscriptRecord {
    pub assembled_transcript_name: Box<str>,
    pub start_chromosome: Box<str>,
    pub start: u32,
    pub end_chromosome: Box<str>,
    pub end: u32,
    pub num_exons: u32,
    pub num_introns: u32
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct DNAVariantRecord {
    pub variant_id: u32,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub strand_1: Box<str>,
    pub operation_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub strand_2: Box<str>,
    pub operation_2: Box<str>,
    pub sequence: Box<str>,
    pub variant_size: u32,
    pub variant_type: Box<str>,
    pub consensus_read_names: Box<str>,
    pub num_consensus_read_names: u32,
    pub read_names: Box<str>,
    pub num_read_names: u32
}

#[derive(Debug, Serialize, Clone)]
pub struct ExonRecord {
    pub assembled_transcript_name: Box<str>,
    pub chromosome: Box<str>,
    pub start: u32,
    pub end: u32,
    pub exon_number: u32,
    pub strand: Box<str>
}

#[derive(Debug, Serialize, Clone)]
pub struct IntronRecord {
    pub assembled_transcript_name: Box<str>,
    pub chromosome: Box<str>,
    pub start: u32,
    pub end: u32,
    pub intron_number: u32,
    pub strand: Box<str>
}

#[derive(Debug, Serialize, Clone)]
pub struct ReferenceTranscriptMatchRecord {
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub num_overlap_bases: u32,
    pub num_transcript_only_bases: u32,
    pub num_reference_transcript_only_bases: u32,
    pub score: f32,
    pub scoring_method: Box<str>
}

#[derive(Debug, Serialize, Clone)]
pub struct ReadFilterStatusRecord {
    pub read_name: Box<str>,
    pub excluded: bool
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TranscriptModelStructureRecord {
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub index: u32,
    pub read_start: u32,
    pub read_end: u32,
    pub sequence: Box<str>,
    #[serde(rename = "type")]
    pub record_type: Box<str>,
    pub kind: Box<str>,
    pub context: Box<str>,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub operation_1: Box<str>,
    pub strand_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub operation_2: Box<str>,
    pub strand_2: Box<str>,
    pub gene_id_1: Box<str>,
    pub transcript_id_1: Box<str>,
    pub exon_id_1: Box<str>,
    pub gene_id_2: Box<str>,
    pub transcript_id_2: Box<str>,
    pub exon_id_2: Box<str>,
    pub skipped: Box<str>
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAVariantRecord {
    pub variant_id: u32,
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub strand_1: Box<str>,
    pub operation_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub strand_2: Box<str>,
    pub operation_2: Box<str>,
    pub variant_size: u32,
    pub variant_type: Box<str>,
    pub sequence: Box<str>,
    pub read_start: u32,
    pub read_end: u32
}
