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

use exacto::caller::prelude as caller;
use exacto::core::prelude as core;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::path::Path;


#[pyfunction]
pub fn identify_somatic_dna_variants(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    control_bam_files: Vec<String>,
    control_bam_bai_files: Vec<String>,
    fasta_file: String,
    output_tsv_file: String,
    regions: Vec<(String, u32, u32)>,
    min_reads: usize,
    min_mapping_quality: u16,
    min_base_quality: u8,
    min_total_depth: u32,
    min_alt_allele_fraction: f32,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    max_slippage_repeat_length: u32,
    num_threads: usize,
    chunk_size: u32,
    max_records: usize,
    expected_variant_allele_fraction: f64,
    expected_mutation_rate: f64,
    expected_sequencing_error: f64,
    expected_slippage_probability: f64,
    max_f1_fraction: f64,
    max_fpr: f64,
    apply_infinite_sites_assumption: bool,
    temp_dir: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let dna_variant_call_set: caller::DNAVariantCallSet = caller::identify_somatic_dna_variants(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        control_bam_files.iter().map(|s| s.as_str()).collect(),
        control_bam_bai_files.iter().map(|s| s.as_str()).collect(),
        fasta_file.as_str(),
        &regions.iter().map(|(chrom, start, end)| (chrom.as_str(), *start, *end)).collect(),
        min_reads,
        min_mapping_quality,
        min_base_quality,
        min_total_depth,
        min_alt_allele_fraction,
        min_size_proportion,
        max_ins_norm_edit_distance,
        max_intrachromosomal_distance_tau,
        max_intrachromosomal_distance,
        max_interchromosomal_distance,
        max_slippage_repeat_length,
        num_threads,
        chunk_size,
        max_records,
        expected_variant_allele_fraction,
        expected_mutation_rate,
        expected_sequencing_error,
        expected_slippage_probability,
        max_f1_fraction,
        max_fpr,
        apply_infinite_sites_assumption,
        temp_dir.as_str()
    );

    match output_type.as_str() {
        "dataframe" => {
            Ok(PyDataFrame(
                caller::dna_variant_records_to_dataframe(
                    caller::build_dna_variant_records(&dna_variant_call_set)
                )
            ))
        }
        "file" => {
            core::write_tsv_file(
                caller::build_dna_variant_records(&dna_variant_call_set),
                Path::new(output_tsv_file.as_str())
            );
            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
