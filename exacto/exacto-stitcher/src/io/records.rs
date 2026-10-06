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
    pub struct StitchedTranscriptRecord {
    pub read_name: Box<str>,
    pub stitched_sequence: Box<str>,
    pub read_length: u32,
    pub stitched_length: u32,

    /// The read bases kept in the stitched sequence, `[retained_read_start, retained_read_end)`
    /// in read coordinates. A 3' roll-back deletes the read bases after `retained_read_end`.
    pub retained_read_start: u32,
    pub retained_read_end: u32,

    /// Five prime
    pub is_five_prime_stitched: bool,
    pub five_prime_reference_gene_id: Box<str>,
    pub five_prime_reference_transcript_id: Box<str>,
    pub five_prime_stitch_start: u32,
    pub five_prime_stitch_end: u32,
    pub five_prime_stitch_length: u32,

    /// Three prime
    pub is_three_prime_stitched: bool,
    pub three_prime_reference_gene_id: Box<str>,
    pub three_prime_reference_transcript_id: Box<str>,
    pub three_prime_stitch_start: u32,
    pub three_prime_stitch_end: u32,
    pub three_prime_stitch_length: u32,

    /// Why the read passed through without a transcript model ("mapping_quality",
    /// "modeling_error"); empty when it was modeled.
    pub unmodeled_reason: Box<str>
}
