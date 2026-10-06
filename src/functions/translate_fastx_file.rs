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

use super::translator_error;


#[pyfunction]
pub fn translate_fastx_file(
    _py: Python,
    fastx_file: String,
    output_fasta_file: String,
    output_tsv_file: String,
    strategy: String,
    start_codons: Vec<String>,
    num_threads: usize
) -> PyResult<()> {
    let translation_strategy: translator::TranslationStrategy = translator::TranslationStrategy::from_str(strategy.as_str())
        .map_err(|_| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Unsupported value for strategy: {:?}.", strategy)))?;
    let start_codons_set: HashSet<&str> = start_codons.iter().map(|s| s.as_str()).collect();
    translator::translate_fastx_file(
        fastx_file.as_str(),
        output_fasta_file.as_str(),
        output_tsv_file.as_str(),
        translation_strategy,
        start_codons_set,
        num_threads
    ).map_err(translator_error)
}
