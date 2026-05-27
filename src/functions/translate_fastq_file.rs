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
extern crate flate2;
extern crate noodles_fastq;
extern crate pyo3;
extern crate tempfile;

use exacto::core::prelude as core;
use exacto::translator::prelude as translator;
use flate2::read::GzDecoder;
use noodles_fastq as fastq;
use pyo3::prelude::*;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::str::FromStr;
use tempfile::Builder as TempFileBuilder;


/// Translate every read in a FASTQ file (plain or gzipped) into primary
/// structures and write the per-PS rows to a TSV file in `temp_dir`.
/// Returns the path to the written TSV.
///
/// Each FASTQ record's name becomes the `Transcript.id` directly, so it
/// surfaces in the output rows as the `transcript_id` column of
/// `PrimaryStructureRecord`.
#[pyfunction]
pub fn translate_fastq_file(
    _py: Python,
    fastq_file: String,
    strategy: String,
    start_codons: Vec<String>,
    num_threads: usize,
    temp_dir: String,
) -> PyResult<String> {
    // Step 1. Open the FASTQ (auto-detect gzip) and collect every record
    // as a (record_name, sequence) pair. The name becomes Transcript.id.
    let gzipped: bool = core::is_gzipped(fastq_file.as_str());
    let file: File = File::open(fastq_file.as_str())
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let reader: Box<dyn Read> = if gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buffered_reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(buffered_reader);

    let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();
    for result in fastq_reader.records() {
        let record = result
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading FASTQ record: {}", e)))?;
        let sequence: String = String::from_utf8(record.sequence().to_vec())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Non-UTF8 sequence: {}", e)))?;
        let id: Box<str> = record.name().to_string().into_boxed_str();
        sequences.push((id, sequence.into_boxed_str()));
    }

    // Step 2. Translate every sequence.
    let translation_strategy: translator::TranslationStrategy =
        translator::TranslationStrategy::from_str(strategy.as_str()).unwrap();
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
    // one PrimaryStructureRecord per row, return its path.
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
        &persisted_file_path,
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    Ok(persisted_file_path.to_string_lossy().to_string())
}
