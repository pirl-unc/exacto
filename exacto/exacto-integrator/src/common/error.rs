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


use thiserror::Error;


/// Why an integration stopped. Each variant names the record that is the cause.
#[derive(Debug, Error)]
pub enum IntegratorError {
    #[error("DNA variant id {variant_id} is held by two calls (origins {origin_1} and {origin_2}); integrate one DNA variant file per origin")]
    DuplicateDNAVariantId { variant_id: u32, origin_1: Box<str>, origin_2: Box<str> },

    #[error("reference transcript {transcript_id} of RNA variant {rna_variant_id} ({assembled_transcript_name}) is not in the gene annotation; use the annotation the RNA variants were called with")]
    UnknownTranscript { transcript_id: Box<str>, rna_variant_id: u32, assembled_transcript_name: Box<str> },

    #[error("a pool of {num_threads} threads could not be built: {reason}")]
    ThreadPool { num_threads: usize, reason: Box<str> }
}
