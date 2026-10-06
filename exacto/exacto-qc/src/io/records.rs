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


#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct TranscriptEndCorrectionRecord {
    pub assembled_transcript_name: Box<str>,
    pub reference_gene_id: Box<str>,
    pub reference_gene_name: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub strand: Box<str>,
    pub end_chromosome: Box<str>,
    pub end_position: u32,
    pub end_context: Box<str>,
    pub end_exon_id: Box<str>,
    pub end_intron_number: u16,
    pub priming_window_sequence: Box<str>,
    pub priming_window_a_count: u32,
    pub priming_window_max_consecutive_a: u32,
    pub upstream_pas_hexamer: Box<str>,
    pub upstream_pas_offset: u32,
    pub num_skipped_reference_bases: u32,
    pub classification: Box<str>,
    pub action: Box<str>,
    pub corrected_read_start: u32,
    pub corrected_read_end: u32,
    pub num_bases_trimmed: u32,
    pub three_prime_complete: bool
}
