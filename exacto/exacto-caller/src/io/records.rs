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
use serde::{Deserialize, Serialize};


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptRecord {
    pub assembled_transcript_name: ReadName,
    pub start_chromosome: ReferenceChromosomeName,
    pub start: ReferencePosition,
    pub end_chromosome: ReferenceChromosomeName,
    pub end: ReferencePosition,
    pub num_exons: u32,
    pub num_splice_junctions: u32,
    pub is_variant: bool
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptExonRecord {
    pub assembled_transcript_name: ReadName,
    pub chromosome: ReferenceChromosomeName,
    pub start: ReferencePosition,
    pub end: ReferencePosition,
    pub exon_number: u32,
    pub strand: Box<str>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptSpliceJunctionRecord {
    pub assembled_transcript_name: ReadName,
    pub chromosome_1: ReferenceChromosomeName,
    pub chromosome_2: ReferenceChromosomeName,
    pub position_1: ReferencePosition,
    pub position_2: ReferencePosition,
    pub strand_1: Box<str>,
    pub strand_2: Box<str>,
    pub splice_junction_number: u32
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptReferenceTranscriptMatchRecord {
    pub assembled_transcript_name: ReadName,
    pub reference_gene_id: ReferenceGeneID,
    pub reference_gene_name: ReferenceGeneName,
    pub reference_transcript_id: ReferenceTranscriptID,
    pub num_splice_junction_matches: u32,
    pub num_overlapping_bases: u32,
    pub num_reference_only_bases: u32,
    pub num_query_only_bases: u32,
    pub num_terminal_query_only_bases: u32,
    pub num_internal_query_only_bases: u32
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptNonsenseMediatedDecayRecord {
    pub assembled_transcript_name: ReadName,
    /// First base of the start codon, as a 0-based read position.
    pub orf_start: ReadPosition,
    
    /// Last base of the stop codon, as a 0-based read position.
    pub orf_end: ReadPosition,
    
    pub nmd_predicted: bool,
    
    /// Exonic bases between the stop codon and the last splice junction. Empty when no junction
    /// lies downstream of the stop codon.
    pub distance_to_last_junction: Option<u32>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptFilterStatusRecord {
    pub read_name: ReadName,
    pub excluded: bool
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptModelAlignmentRecord {
    pub assembled_transcript_name: ReadName,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub index: u32,
    pub read_start: ReadPosition,
    pub read_end: ReadPosition,
    pub sequence: Box<str>,
    #[serde(rename = "type")]
    pub record_type: Box<str>,
    pub kind: Box<str>,
    pub context: Box<str>,
    pub chromosome_1: ReferenceChromosomeName,
    pub position_1: ReferencePosition,
    pub operation_1: Box<str>,
    pub strand_1: Box<str>,
    pub chromosome_2: ReferenceChromosomeName,
    pub position_2: ReferencePosition,
    pub operation_2: Box<str>,
    pub strand_2: Box<str>,
    pub reference_gene_id_1: ReferenceGeneID,
    pub reference_transcript_id_1: ReferenceTranscriptID,
    pub reference_exon_id_1: ReferenceExonID,
    pub reference_gene_id_2: ReferenceGeneID,
    pub reference_transcript_id_2: ReferenceTranscriptID,
    pub reference_exon_id_2: ReferenceExonID,
    pub skipped: Box<str>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptVariantRecord {
    pub variant_id: u32,
    pub assembled_transcript_name: ReadName,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub chromosome_1: ReferenceChromosomeName,
    pub position_1: ReferencePosition,
    pub strand_1: Box<str>,
    pub operation_1: Box<str>,
    pub chromosome_2: ReferenceChromosomeName,
    pub position_2: ReferencePosition,
    pub strand_2: Box<str>,
    pub operation_2: Box<str>,
    
    /// Empty where a size has no meaning: a translocation, a fusion or a circular RNA.
    pub variant_size: Option<u32>,
    
    pub variant_type: Box<str>,
    
    /// Inserted or substituted bases in forward (reference) orientation, whatever the strand.
    pub sequence: Box<str>,
    
    pub read_start: ReadPosition,
    pub read_end: ReadPosition,
    
    /// Origins of the DNA variants the call matches (`germline`, `somatic`), `;`-joined in the
    /// order of the DNA variant tables; empty for no match or no DNA variants given.
    #[serde(default)]
    pub origin: Box<str>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct DNAVariantRecord {
    pub variant_id: u32,
    pub origin: Box<str>,
    pub chromosome_1: ReferenceChromosomeName,
    pub position_1: ReferencePosition,
    pub strand_1: Box<str>,
    pub operation_1: Box<str>,
    pub chromosome_2: ReferenceChromosomeName,
    pub position_2: ReferencePosition,
    pub strand_2: Box<str>,
    pub operation_2: Box<str>,
    
    /// Inserted or substituted bases in forward (reference) orientation, whatever the strand.
    pub sequence: Box<str>,
    
    /// Empty where a size has no meaning: a translocation.
    #[serde(default)]
    pub variant_size: Option<i32>,
    
    #[serde(default)]
    pub variant_type: Box<str>,
    #[serde(default)]
    pub consensus_read_names: Box<str>,
    #[serde(default)]
    pub num_consensus_read_names: u32,
    #[serde(default)]
    pub read_names: Box<str>,
    #[serde(default)]
    pub num_read_names: ReadSupport
}


#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterRecord {
    pub cluster_id: ClusterID,
    pub read_name: ReadName
}
