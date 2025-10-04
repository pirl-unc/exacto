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
extern crate polars;
extern crate pyo3;
extern crate pyo3_polars;

use exacto::annotator::prelude as annotator;
use exacto::core::prelude as core;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use std::collections::HashSet;


#[pyfunction]
pub fn annotate_variant_calls(
    py: Python,
    tsv_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    gene_types: Vec<String>,
    gene_levels: Vec<u8>,
    transcript_types: Vec<String>,
    transcript_levels: Vec<u8>,
    output_tsv_file: String,
    num_threads: usize,
    temp_dir: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let gene_annotator = if reference_gene_annotation_source.as_str() == "gencode" {
        let gene_types_: Option<HashSet<&str>> = (!gene_types.is_empty()).then(|| gene_types.iter().map(String::as_str).collect());
        let gene_levels_: Option<HashSet<u8>> = (!gene_levels.is_empty()).then(|| gene_levels.iter().copied().collect());
        let transcript_types_: Option<HashSet<&str>> = (!transcript_types.is_empty()).then(|| transcript_types.iter().map(String::as_str).collect());
        let transcript_levels_: Option<HashSet<u8>> = (!transcript_levels.is_empty()).then(|| transcript_levels.iter().copied().collect());
        core::Gencode::new(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str(),
            gene_types_,
            gene_levels_,
            transcript_types_,
            transcript_levels_
        )
    } else {
        panic!("Unsupported annotation source: {}", reference_gene_annotation_source);
    };
    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();
    let variant_call_annotation_set: annotator::VariantCallAnnotationSet = annotator::annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        num_threads
    );
    match output_type.as_str() {
        "dataframe" => {
            Ok(PyDataFrame(variant_call_annotation_set.to_dataframe()))
        }
        "file" => {
            variant_call_annotation_set.to_tsv_file(
                output_tsv_file.as_str()
            );
            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        },
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
