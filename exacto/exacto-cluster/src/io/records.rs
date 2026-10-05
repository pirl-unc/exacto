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


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterIDRecord {
    pub cluster_id: usize,
    pub read_name: Box<str>,
    
    /// The read was assigned to another cluster and fits this one as well. A table written
    /// before the column existed reads as false.
    #[serde(default)]
    pub is_shared: bool
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterSummaryRecord {
    pub cluster_id: usize,
    pub num_reads: usize
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterReferenceTranscriptRecord {
    pub cluster_id: usize,
    pub reference_gene_id: Box<str>,
    pub reference_transcript_id: Box<str>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterSpliceJunctionRecord {
    pub cluster_id: usize,
    pub chromosome_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_1: u32,
    pub position_2: u32,
    pub strand_1: Box<str>,
    pub strand_2: Box<str>
}
#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterVariantRecord {
    pub cluster_id: usize,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub strand_1: Box<str>,
    pub operation_type_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub strand_2: Box<str>,
    pub operation_type_2: Box<str>,
    pub sequence: Box<str>,
    pub variant_type: Box<str>
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterTemplateSwitchRecord {
    pub cluster_id: usize,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub strand_1: Box<str>,
    pub operation_type_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub strand_2: Box<str>,
    pub operation_type_2: Box<str>,
    pub sequence: Box<str>,
    pub variant_type: Box<str>,
    pub homology_left: Option<u32>,
    pub homology_right: Option<u32>,
    pub canonical_splice: Option<bool>,
    pub breakpoint_dispersion: Option<u32>,
    pub foldback: bool,
    pub partner_count: Option<u32>,
    pub num_reads: u32,
    pub flagged: bool,
    pub suppressed: bool
}


#[derive(Debug, Default, Serialize, Deserialize, Clone, PartialEq)]
pub struct RNAReadClusterFailedVariantRecord {
    pub cluster_id: Option<usize>,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub strand_1: Box<str>,
    pub operation_type_1: Box<str>,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub strand_2: Box<str>,
    pub operation_type_2: Box<str>,
    pub sequence: Box<str>,
    pub variant_type: Box<str>,
    pub failure_reason: Box<str>,
    pub alt_count: u32,
    pub reference_count: u32,
    pub coverage: Option<u32>,
    pub p_value: Option<f64>,
    pub p_value_cutoff: Option<f64>,
    pub expected_error_rate: Option<f64>,
    pub slippage_repeat_length: Option<u32>,
    #[serde(default)]
    pub junction_homology: Option<u32>
}
