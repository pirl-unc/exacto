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


extern crate chrono;
extern crate env_logger;
extern crate exitcode;
extern crate log;
extern crate pyo3;
extern crate serde_json;
use chrono::Local;
use env_logger::{Builder, Env};
use log::{info, LevelFilter};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::io::Write;
mod clustering;
mod constants;
mod defaults;
mod variant_identification_dna;
mod variant_identification_rna;
mod utilities;
mod variant_call;
mod variant_record;
use variant_identification_dna::identify_dna_variants;
use variant_identification_rna::identify_rna_variants;
use variant_call::VariantCall;
use variant_record::VariantRecord;


/// This function identifies DNA variants in a long-read WGS BAM file.
///
/// # Arguments
/// * `bam_file`                    -   BAM file.
/// * `min_reads`                   -   Minimum number of reads.
/// * `min_mapping_quality`         -   Minimum mapping quality.
/// * `num_threads`                 -   number of threads.
/// * `chromosomes`                 -   Chromosomes to call.
///
/// # Returns
/// * `nearby_variants_map`         -   HashMap where key is Variant.id and
///                                     value is a vector of GenomicRange.id
#[pyfunction]
fn call_dna_variants(
    py: Python,
    bam_file: String,
    sample_id: String,
    min_reads: usize,
    min_mapping_quality: usize,
    min_ins_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    min_del_size_proportion: f32,
    max_bnd_distance: isize,
    clustering_grid_size: isize,
    num_threads: usize,
    chromosomes: Vec<String>) -> PyResult<String> {
    let variant_calls: Vec<VariantCall> = identify_dna_variants(
        &bam_file,
        &sample_id,
        min_reads,
        min_mapping_quality,
        num_threads,
        min_ins_size_proportion,
        max_ins_norm_edit_distance,
        min_del_size_proportion,
        max_bnd_distance,
        clustering_grid_size,
        chromosomes.iter().map(|s| s.as_str()).collect()
    );
    let serialized = serde_json::to_string(&variant_calls).expect("Serialization of vector of VariantCall object failed");
    Ok(serialized)
}


/// This function identifies RNA variants in a long-read RNA-seq BAM file.
///
/// # Arguments
/// * `bam_file`                    -   BAM file.
/// * `min_reads`                   -   Minimum number of reads.
/// * `min_mapping_quality`         -   Minimum mapping quality.
/// * `num_threads`                 -   number of threads.
/// * `chromosomes`                 -   Chromosomes to call.
///
/// # Returns
/// * `nearby_variants_map`         -   HashMap where key is Variant.id and
///                                     value is a vector of GenomicRange.id
#[pyfunction]
fn call_rna_variants(
    py: Python,
    bam_file: String,
    sample_id: String,
    min_reads: usize,
    min_mapping_quality: usize,
    min_ins_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    min_del_size_proportion: f32,
    max_bnd_distance: isize,
    clustering_grid_size: isize,
    num_threads: usize,
    chromosomes: Vec<String>) -> PyResult<String> {
    let variant_calls: Vec<VariantCall> = identify_rna_variants(
        &bam_file,
        &sample_id,
        min_reads,
        min_mapping_quality,
        num_threads,
        min_ins_size_proportion,
        max_ins_norm_edit_distance,
        min_del_size_proportion,
        max_bnd_distance,
        clustering_grid_size,
        chromosomes.iter().map(|s| s.as_str()).collect()
    );
    let serialized = serde_json::to_string(&variant_calls).expect("Serialization of vector of VariantCall object failed");
    Ok(serialized)
}


#[pymodule]
fn exactolibrs(_py: Python, m: &PyModule) -> PyResult<()> {
    // Initialize the logger
    Builder::from_env(Env::default().default_filter_or("info")).format(|buf, record| {
        let now = Local::now();
        writeln!(
            buf,
            "{} {} [{:>50}] {}",
            now.format("%Y-%m-%d %H:%M:%S"),
            record.level(),
            record.target(),
            record.args()
        )
    }).init();

    m.add_function(wrap_pyfunction!(call_dna_variants, m)?);
    m.add_function(wrap_pyfunction!(call_rna_variants, m)?);
    Ok(())
}
