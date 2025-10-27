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

use exacto::caller::prelude as caller;
use exacto::core::prelude as core;
use exacto::translator::prelude as translator;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use std::collections::HashSet;
use std::str::FromStr;


#[pyfunction]
pub fn translate_structures(
    py: Python,
    transcript_structures_tsv_file: String,
    rna_variant_calls_tsv_file: String,
    integrated_variants_tsv_file: String,
    strategy: String,
    output_tsv_file: String,
    num_threads: usize,
    output_type: String
) -> PyResult<PyDataFrame> {
    let df_transcript_structures: DataFrame = core::read_tsv_file(transcript_structures_tsv_file.as_str());
    let rna_variant_call_set: caller::RNAVariantCallSet = caller::RNAVariantCallSet::read_tsv_file(rna_variant_calls_tsv_file.as_str());
    let df_integrated_variants: DataFrame = core::read_tsv_file(integrated_variants_tsv_file.as_str());
    let strategy_: translator::TranslationStrategy = translator::TranslationStrategy::from_str(strategy.as_str()).unwrap();

    let primary_structure_set: translator::PrimaryStructureSet = translator::translate_transcript_structures(
        &df_transcript_structures,
        &rna_variant_call_set,
        &df_integrated_variants,
        strategy_,
        num_threads
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_primary_structures: DataFrame = primary_structure_set.to_dataframe();
            Ok((PyDataFrame(df_primary_structures)))
        },
        "file" => {
            primary_structure_set.to_tsv_file(output_tsv_file.as_str());
            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap())))
        },
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
