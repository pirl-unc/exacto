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
use std::collections::HashSet;


#[pyfunction]
pub fn identify_rna_variants(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    reference_genome_fasta_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    gene_types: Vec<String>,
    gene_levels: Vec<u8>,
    transcript_types: Vec<String>,
    transcript_levels: Vec<u8>,
    output_dir: String,
    output_prefix: String,
    reference_transcript_scoring_method: String,
    reference_transcript_selection_strategy: String,
    reference_transcript_top_k: usize,
    reference_transcript_threshold: f32,
    min_mapping_quality: u32,
    min_average_base_quality: u8,
    num_threads: usize,
    temp_dir: String,
    output_type: String
) -> PyResult<(PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame)> {
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
    let scoring_method: caller::ReferenceTranscriptScoringMethod = reference_transcript_scoring_method.as_str().parse().unwrap();
    let selection_strategy: caller::ReferenceTranscriptSelectionStrategy = reference_transcript_selection_strategy.as_str().parse().unwrap();
    let transcript_model_set: caller::TranscriptModelSet = caller::identify_variant_transcripts(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        reference_genome_fasta_file.as_str(),
        &gene_annotator,
        scoring_method,
        selection_strategy,
        reference_transcript_top_k,
        reference_transcript_threshold,
        min_mapping_quality,
        min_average_base_quality,
        num_threads
    );
    match output_type.as_str() {
        "dataframe" => {
            let df_exons = transcript_model_set.get_exons_dataframe();
            let df_read_filter_status: DataFrame = transcript_model_set.get_read_filter_status_dataframe();
            let df_read_names: DataFrame = transcript_model_set.get_read_names_dataframe();
            let df_matched_reference_transcripts: DataFrame = transcript_model_set.get_reference_transcript_matches_dataframe();
            let df_introns: DataFrame = transcript_model_set.get_introns_dataframe();
            let df_transcripts: DataFrame = transcript_model_set.get_transcripts_dataframe();
            let df_variant_calls: DataFrame = transcript_model_set.get_variant_calls_dataframe(num_threads);
            Ok((PyDataFrame(df_exons),
                PyDataFrame(df_read_filter_status),
                PyDataFrame(df_read_names),
                PyDataFrame(df_matched_reference_transcripts),
                PyDataFrame(df_introns),
                PyDataFrame(df_transcripts),
                PyDataFrame(df_variant_calls)))
        }
        "file" => {
            transcript_model_set.to_tsv_files(
                output_dir.as_str(),
                output_prefix.as_str()
            );
            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap())))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
