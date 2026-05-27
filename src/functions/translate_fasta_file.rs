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
extern crate tempfile;

use exacto::core::prelude as core;
use exacto::translator::prelude as translator;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::env;
use std::path::Path;
use std::str::FromStr;
use tempfile::Builder as TempFileBuilder;


#[pyfunction]
pub fn translate_fasta_file(
    _py: Python,
    fasta_file: String,
    strategy: String,
    start_codons: Vec<String>,
    num_threads: usize,
    temp_dir: String
) -> PyResult<String> {
    // Step 1. Read every FASTA record into (sequence_id, sequence) pairs.
    // The FASTA sequence name becomes Transcript.id directly.
    let sequence_ids: Vec<(Box<str>, u32)> = core::get_fasta_sequence_ids(fasta_file.as_str());
    let sequences: Vec<(Box<str>, Box<str>)> = sequence_ids
        .into_iter()
        .map(|(sequence_id, length)| {
            let sequence: Box<str> = core::get_fasta_sequence(&*sequence_id, 1, length, fasta_file.as_str());
            (sequence_id, sequence)
        })
        .collect();

    // Step 2. Translate every sequence.
    let translation_strategy: translator::TranslationStrategy = translator::TranslationStrategy::from_str(strategy.as_str()).unwrap();
    let start_codons_set: HashSet<&str> = start_codons.iter().map(|s| s.as_str()).collect();
    let transcript_set: translator::TranscriptSet = translator::translate_sequences(
        sequences,
        translation_strategy,
        start_codons_set,
        num_threads,
    );

    // Step 3. Resolve the output directory.
    let dir: String = if temp_dir.is_empty() {
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Directory does not exist: {}", dir),
        ));
    }

    // Step 4. Create a uniquely-named .tsv file in the directory, write
    // one PrimaryStructureLevelRecord per row, return its path.
    let temp_file = TempFileBuilder::new()
        .suffix(".tsv")
        .tempfile_in(dir_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let persisted_file_path = temp_file
        .into_temp_path()
        .keep()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    core::write_tsv_file(
        translator::build_primary_structure_records(&transcript_set),
        &persisted_file_path
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    Ok(persisted_file_path.to_string_lossy().to_string())
}
