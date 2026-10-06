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


use csv::ReaderBuilder;

use crate::prelude::*;


pub fn load_consensus_sequence_records(tsv_file: &str) -> Vec<ConsensusSequenceRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .unwrap_or_else(|error| panic!("Failed to open the RNA consensus TSV file {}: {}", tsv_file, error));
    reader
        .deserialize()
        .map(|result| result.unwrap_or_else(|error| panic!("Failed to read a row of the RNA consensus TSV file {}: {}", tsv_file, error)))
        .collect()
}


#[cfg(test)]
#[path = "../tests/io/loaders.rs"]
mod tests;