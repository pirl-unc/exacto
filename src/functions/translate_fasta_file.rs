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
extern crate noodles_fasta;
extern crate polars;
extern crate pyo3;
extern crate pyo3_polars;
extern crate tempfile;

use exacto::core::prelude as core;
use exacto::translator::prelude as translator;
use flate2::read::GzDecoder;
use noodles_fasta as fasta;
use noodles_fasta::{fai, Record};
use noodles_fasta::io::indexed_reader::Builder;
use noodles_fasta::io::Reader as FastaReader;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufReader,BufWriter,Cursor,Read,Write};
use std::path::Path;
use std::str::FromStr;
use tempfile::NamedTempFile;


#[pyfunction]
pub fn translate_fasta_file(
    py: Python,
    fasta_file: String,
    strategy: String,
    num_threads: usize,
    temp_dir: String
) -> PyResult<String> {
    // Step 1. Read the FASTA file
    let mut rnas: Vec<translator::RNA> = Vec::new();
    let sequence_ids: Vec<(Box<str>, u32)> = core::get_fasta_sequence_ids(fasta_file.as_str());
    for (sequence_id, length) in sequence_ids.iter() {
        let sequence: Box<str> = core::get_fasta_sequence(&**sequence_id, 1, *length, fasta_file.as_str());
        let rna = translator::RNA::new(sequence_id.clone(), sequence);
        rnas.push(rna);
    }

    // Step 2. Translate the RNA sequences
    let mut translation_set: translator::TranslationSet = translator::translate_rnas(
        rnas,
        num_threads
    );

    // Step 3. Write the translated peptides set to a temporary file
    let dir: String = if temp_dir.is_empty() {
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }
    let temp_file = NamedTempFile::new_in(dir_path).unwrap();
    let persisted_file_path = temp_file.into_temp_path().keep().unwrap();
    let file = File::create(&persisted_file_path).unwrap();
    let mut writer = IpcWriter::new(file);
    let strategy: translator::TranslationStrategy = translator::TranslationStrategy::from_str(strategy.as_str()).unwrap();
    writer.finish(&mut translation_set.to_dataframe(strategy)).unwrap();
    drop(writer);
    translation_set.translations.clear();
    translation_set.translations.shrink_to_fit();
    drop(translation_set);
    Ok(persisted_file_path.to_str().unwrap().to_string())
}
