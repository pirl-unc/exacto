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


extern crate bam;
mod variant_calling;
mod utils;

use pyo3::prelude::*;
use std::str;
use lazy_static::lazy_static;
use regex::Regex;
use polars::prelude::*;
use polars::df;
use rustc_hash::FxHashMap;
use variant_calling::*;
use utils::*;


#[pyclass]
struct VariantCallset {
    #[pyo3(get)]
    read_ids: Vec<String>,
    #[pyo3(get)]
    chromosomes: Vec<String>,
    #[pyo3(get)]
    positions: Vec<i32>,
    #[pyo3(get)]
    variant_types: Vec<String>,
    #[pyo3(get)]
    reference_alleles: Vec<String>,
    #[pyo3(get)]
    alternate_alleles: Vec<String>,
    #[pyo3(get)]
    sequences: Vec<String>,
    #[pyo3(get)]
    variant_sizes: Vec<i32>
}

#[pymethods]
impl VariantCallset {
    #[new]
    pub fn new(
        read_ids: Vec<String>,
        chromosomes: Vec<String>,
        positions: Vec<i32>,
        variant_types: Vec<String>,
        reference_alleles: Vec<String>,
        alternate_alleles: Vec<String>,
        sequences: Vec<String>,
        variant_sizes: Vec<i32>) -> Self {
        VariantCallset {
            read_ids: read_ids.clone(),
            chromosomes: chromosomes.clone(),
            positions: positions.clone(),
            variant_types: variant_types.clone(),
            reference_alleles: reference_alleles.clone(),
            alternate_alleles: alternate_alleles.clone(),
            sequences: sequences.clone(),
            variant_sizes: variant_sizes.clone()
        }
    }
}

/// Identifies RNA variants in a BAM file.
#[pyfunction]
fn identify_rna_variants(bam_file: String, num_threads: u16) -> PyResult<VariantCallset> {
    let mut df_cs_tags = DataFrame::default();
    let chromosomes: Vec<String> = get_chromosome_names(&bam_file);
    let reader = bam::BamReader::from_path(bam_file, num_threads).unwrap();
    for record in reader {
        let record = record.unwrap();
        let record_id: &str = str::from_utf8(&record.name()).unwrap();
        let ref_id: usize = record.ref_id().try_into().unwrap();
        let chromosome: &str = &chromosomes[ref_id];
        let flag = record.flag();

        // Check if the read is mapped
        if flag.is_mapped() == false {
            println!("Unmapped read ID: {}", record_id);
            continue;
        }

        let start_pos: i32 = record.start(); // originally 0-based
        let mapping_quality: u8 = record.mapq();
        let sequence_string: String = record.sequence().to_vec_acgtn_only().iter().map(|i| (*i as char)).collect();
        let sequence_vec: Vec<char> = record.sequence().to_vec_acgtn_only().iter().map(|i| (*i as char)).collect();
        let cigar = record.cigar();
        let mut cs_tag: &str = "";
        match record.tags().get(b"cs") {
            Some(bam::record::tags::TagValue::String(tag_value, bam::record::tags::StringType::String)) => {
                cs_tag = str::from_utf8(tag_value).unwrap();
            },
            Some(bam::record::tags::TagValue::Char(value)) => println!("Char = {}", value),
            _ => panic!("Unexpected type"),
        }

        let df_curr_cs_tags: DataFrame = identify_rna_variants_in_cs_tag(record_id, chromosome, start_pos, cs_tag).unwrap();
        df_cs_tags = df_cs_tags.vstack(&df_curr_cs_tags).unwrap();

        // println!("Read ID: {}", record_id);
        // println!("Mapping quality: {}", mapping_quality);
        // println!("Sequence (string): {}", sequence_string);
        // println!("Sequence (vector): {:?}", sequence_vec);
        // println!("First character of sequence: {:?}", sequence_vec[0]);
        // println!("First CIGAR: {:?}", cigar.at(0));
        // println!("CIGAR length: {}", cigar.len());
        // println!("CS tag: {}", cs_tag);
    }

    // Sort the DataFrame
    df_cs_tags = df_cs_tags.sort(&["position"], vec![false, true, false, false, false, false, false, false]).unwrap();

    // Convert the DataFrame to a wrapper
    let read_ids: Vec<String> = copy_string_series_as_vector(&df_cs_tags["read_id"]);
    let chromosomes: Vec<String> = copy_string_series_as_vector(&df_cs_tags["chromosome"]);
    let positions: Vec<i32> = copy_i32_series_as_vector(&df_cs_tags["position"]);
    let variant_types: Vec<String> = copy_string_series_as_vector(&df_cs_tags["variant_type"]);
    let reference_alleles: Vec<String> = copy_string_series_as_vector(&df_cs_tags["reference_allele"]);
    let alternate_alleles: Vec<String> = copy_string_series_as_vector(&df_cs_tags["alternate_allele"]);
    let sequences: Vec<String> = copy_string_series_as_vector(&df_cs_tags["sequence"]);
    let variant_sizes: Vec<i32> = copy_i32_series_as_vector(&df_cs_tags["variant_size"]);
    let variant_callset = VariantCallset {
        read_ids: read_ids,
        chromosomes: chromosomes,
        positions: positions,
        variant_types: variant_types,
        reference_alleles: reference_alleles,
        alternate_alleles: alternate_alleles,
        sequences: sequences,
        variant_sizes: variant_sizes
    };
    Ok(variant_callset)
}

/// A Python module implemented in Rust.
#[pymodule]
fn exactors(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<VariantCallset>()?;
    m.add_function(wrap_pyfunction!(identify_rna_variants, m)?)?;
    Ok(())
}