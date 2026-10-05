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


extern crate exacto;
extern crate pyo3;

use exacto::translator::prelude as translator;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::str::FromStr;

use super::translator_error;


/// Translate a single RNA sequence and return one `(proteoform_sequence,
/// orf_start, orf_end)` tuple per proteoform produced.
///
/// `strategy` selects how many ORFs come back (`LongestORF` → at most one;
/// `AllORFs` → all viable ORFs). The lone synthetic transcript is built
/// with id `"sequence"` so this function stays as light as possible.
#[pyfunction]
pub fn translate_sequence(
    _py: Python,
    rna_sequence: String,
    strategy: String,
    start_codons: Vec<String>,
) -> PyResult<Vec<(String, u32, u32)>> {
    let translation_strategy: translator::TranslationStrategy =
        translator::TranslationStrategy::from_str(strategy.as_str())
            .map_err(|_| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Unsupported value for strategy: {:?}.", strategy)))?;
    let start_codons_set: HashSet<&str> = start_codons.iter().map(|s| s.as_str()).collect();

    let transcript_set: translator::AssembledTranscriptSet = translator::translate_sequences(
        vec![("sequence".to_string().into_boxed_str(), rna_sequence.into_boxed_str())],
        translation_strategy,
        &start_codons_set,
        1,
    ).map_err(translator_error)?;

    // Walk the one Transcript's proteoform
    let mut proteoforms: Vec<(String, u32, u32)> = Vec::new();
    for transcript in transcript_set.iter() {
        for proteoform in transcript.proteoforms.iter() {
            proteoforms.push((proteoform.get_sequence().to_string(), proteoform.get_orf_start(), proteoform.get_orf_end()));
        }
    }
    Ok(proteoforms)
}
