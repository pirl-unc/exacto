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

use crate::functions::io_error;


#[pyfunction]
pub fn identify_germline_dna_variants(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    fasta_file: String,
    output_tsv_file: String,
    regions: Vec<(String, u32, u32)>,
    min_reads: usize,
    min_mapping_quality: u16,
    min_base_quality: u8,
    min_total_depth: usize,
    min_alt_allele_fraction: f64,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_intrachromosomal_distance: u32,
    max_slippage_repeat_length: u32,
    num_threads: usize,
    chunk_size: u32,
    max_records: u32,
    expected_variant_allele_fraction: f64,
    expected_mutation_rate: f64,
    expected_sequencing_error: f64,
    expected_slippage_probability: f64,
    max_fpr: f64,
    min_terminal_softclip_length: u32,
    bkpt_rescue: bool,
    bkpt_rescue_min_ins_len: u32,
    bkpt_rescue_max_ins_len: u32,
    bkpt_rescue_search_distance: u32,
    bkpt_rescue_realignment_gap_open_score: i32,
    bkpt_rescue_realignment_gap_extend_score: i32,
    bkpt_rescue_realignment_k: u32,
    bkpt_rescue_realignment_band_width: u32,
    bkpt_rescue_realignment_min_score_fraction: f64,
    bkpt_rescue_realignment_min_query_coverage: f64,
    bkpt_rescue_realignment_min_span_proportion: f64,
    poa_match_score: i32,
    poa_mismatch_score: i32,
    poa_gap_open_score: i32,
    poa_gap_extend_score: i32,
    read_depth_max_merge_distance: u32,
    min_homopolymer_len: u32,
    min_dinucleotide_context_len: u32,
    temp_dir: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let _ = py;
    // Every option field is a parameter, so the Python defaults are the only defaults.
    let options: caller::IdentifyGermlineDNAVariantsOptions = caller::IdentifyGermlineDNAVariantsOptions {
        error_model: caller::SequencingErrorModel {
            sequencing_error: expected_sequencing_error,
            slippage_prob: expected_slippage_probability
        },
        prior: caller::DNAVariantPriorModel {
            allele_fraction: expected_variant_allele_fraction,
            mutation_rate: expected_mutation_rate
        },
        calling: caller::DNAVariantCallingOptions {
            chunk_size: chunk_size,
            max_records: max_records,
            min_mapping_quality: min_mapping_quality,
            min_base_quality: min_base_quality,
            max_ins_norm_edit_distance: max_ins_norm_edit_distance,
            min_size_proportion: min_size_proportion,
            max_clustering_distance: max_intrachromosomal_distance,
            min_terminal_soft_clip_ins_len: min_terminal_softclip_length,
            bkpt_rescue: bkpt_rescue,
            bkpt_rescue_min_ins_len: bkpt_rescue_min_ins_len,
            bkpt_rescue_max_ins_len: bkpt_rescue_max_ins_len,
            bkpt_rescue_search_distance: bkpt_rescue_search_distance,
            bkpt_rescue_realignment_gap_open_score: bkpt_rescue_realignment_gap_open_score,
            bkpt_rescue_realignment_gap_extend_score: bkpt_rescue_realignment_gap_extend_score,
            bkpt_rescue_realignment_k: bkpt_rescue_realignment_k,
            bkpt_rescue_realignment_band_width: bkpt_rescue_realignment_band_width,
            bkpt_rescue_realignment_min_score_fraction: bkpt_rescue_realignment_min_score_fraction,
            bkpt_rescue_realignment_min_query_coverage: bkpt_rescue_realignment_min_query_coverage,
            bkpt_rescue_realignment_min_span_proportion: bkpt_rescue_realignment_min_span_proportion,
            poa_match_score: poa_match_score,
            poa_mismatch_score: poa_mismatch_score,
            poa_gap_open_score: poa_gap_open_score,
            poa_gap_extend_score: poa_gap_extend_score,
            read_depth_max_merge_distance: read_depth_max_merge_distance
        },
        filtering: caller::DNAVariantFilteringOptions {
            max_fpr: max_fpr,
            min_reads: min_reads,
            min_total_depth: min_total_depth,
            min_alt_allele_fraction: min_alt_allele_fraction,
            max_slippage_repeat_len: max_slippage_repeat_length,
            min_homopolymer_len: min_homopolymer_len,
            min_dinucleotide_context_len: min_dinucleotide_context_len
        }
    };

    let dna_variant_call_set: caller::DNAVariantCallSet = caller::identify_germline_dna_variants(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        fasta_file.as_str(),
        &regions.iter().map(|(chrom, start, end)| (chrom.as_str(), *start, *end)).collect(),
        &options,
        num_threads,
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
            core::write_tsv_table(
                caller::build_dna_variant_records(&dna_variant_call_set),
                Path::new(output_tsv_file.as_str())
            ).map_err(io_error)?;
            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}