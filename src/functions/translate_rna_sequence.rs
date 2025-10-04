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
extern crate polars;
extern crate pyo3;
extern crate pyo3_polars;
extern crate tempfile;

use exacto::core::prelude as core;
use exacto::translator::prelude as translator;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use std::collections::HashSet;
use std::str::FromStr;


#[pyfunction]
pub fn translate_rna_sequence(
    py: Python,
    rna_sequence: String,
    strategy: String
) -> PyResult<Vec<(String, usize, usize)>> {
    let rnas = vec![translator::RNA::new("sequence".into(), rna_sequence.into_boxed_str())];
    let translation_set: translator::TranslationSet = translator::translate_rnas(rnas, 1);
    let strategy_: translator::TranslationStrategy = translator::TranslationStrategy::from_str(strategy.as_str()).unwrap();
    let mut translations: Vec<(String, usize, usize)> = Vec::new();
    for translation in translation_set.get_peptides_by_strategy(strategy_).iter() {
        for peptide in translation.peptides.iter() {
            translations.push((peptide.sequence.to_string(), peptide.orf_start, peptide.orf_end));
        }
    }
    Ok(translations)
}
