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
use exacto::caller::prelude as caller;
use exacto::core::prelude as core;
use exacto::integrator::prelude as integrator;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;


/// This function integrates DNA and RNA variants.
#[pyfunction]
pub fn integrate_dna_rna_variants(
    py: Python,
    dna_variant_call_annotation_set_tsv_file: String,
    rna_variant_call_set_tsv_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    output_tsv_file: String,
    max_exon_offset: u32,
    max_transcript_boundary_offset: usize,
    max_intergenic_distance: usize,
    num_threads: usize,
    output_type: String
) -> PyResult<PyDataFrame> {
    let dna_variant_call_annotation_set: annotator::VariantCallAnnotationSet = annotator::VariantCallAnnotationSet::read_tsv_file(dna_variant_call_annotation_set_tsv_file.as_str());
    let rna_variant_call_set: caller::RNAVariantCallSet = caller::RNAVariantCallSet::read_tsv_file(rna_variant_call_set_tsv_file.as_str());

    let gene_annotator = if reference_gene_annotation_source.as_str() == "gencode" {
        core::Gencode::new_with_defaults(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str()
        )
    } else {
        panic!("Unsupported annotation source: {}", reference_gene_annotation_source);
    };

    let integrated_variant_set: integrator::IntegratedVariantSet = integrator::integrate_dna_rna_variants(
        &dna_variant_call_annotation_set,
        &rna_variant_call_set,
        &gene_annotator,
        max_exon_offset,
        max_transcript_boundary_offset,
        max_intergenic_distance,
        num_threads
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_integrations: DataFrame = integrated_variant_set.to_dataframe();
            Ok((PyDataFrame(df_integrations)))
        },
        "file" => {
            integrated_variant_set.to_tsv_file(
                output_tsv_file.as_str()
            );
            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap())))
        },
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
