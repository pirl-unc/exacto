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

use crate::io::records::*;


/// Loads DNA variant records from a tab-separated file.
///
/// Accepts either the full exacto somatic-variants schema (16 columns) or a
/// bare-bones variant-grammar schema providing only the first 11 columns
/// (`variant_id`, `chromosome_1`, `position_1`, `strand_1`, `operation_1`,
/// `chromosome_2`, `position_2`, `strand_2`, `operation_2`, `sequence`, `origin`). csv maps
/// columns by header name, and the trailing `DNAVariantRecord` fields are
/// `#[serde(default)]`, so a bare-bones file's missing columns default to empty
/// strings / `0`. Column order is irrelevant; extra columns are ignored.
pub fn load_dna_variant_records(tsv_file: &str) -> Vec<DNAVariantRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .expect("Failed to open TSV file");
    reader
        .deserialize()
        .map(|result| result.expect("Failed to deserialize DNAVariantRecord row"))
        .collect()
}


pub fn load_assembled_transcript_variant_records(tsv_file: &str) -> Vec<AssembledTranscriptVariantRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .expect("Failed to open TSV file");
    reader
        .deserialize()
        .map(|result| result.expect("Failed to deserialize RNAVariantRecord row"))
        .collect()
}


pub fn load_assembled_transcript_model_alignment_records(tsv_file: &str) -> Vec<AssembledTranscriptModelAlignmentRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .expect("Failed to open TSV file");
    reader
        .deserialize()
        .map(|result| result.expect("Failed to deserialize TranscriptModelAlignmentRecord row"))
        .collect()
}


pub fn load_assembled_transcript_reference_transcript_match_records(tsv_file: &str) -> Vec<AssembledTranscriptReferenceTranscriptMatchRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .expect("Failed to open TSV file");
    reader
        .deserialize()
        .map(|result| result.expect("Failed to deserialize AssembledTranscriptReferenceTranscriptMatchRecord row"))
        .collect()
}


#[cfg(test)]
#[path = "../tests/io/loaders.rs"]
mod tests;