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
use std::path::{Path, PathBuf};


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
    output_dir: String,
    output_prefix: String,
    reference_transcript_scoring_method: String,
    reference_transcript_selection_strategy: String,
    reference_transcript_top_k: usize,
    reference_transcript_threshold: f32,
    min_mapping_quality: u16,
    min_average_base_quality: u8,
    num_threads: usize,
    temp_dir: String,
    output_type: String,
    chunk_size: usize
) -> PyResult<(PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame)> {
    let gene_annotator = if reference_gene_annotation_source.as_str() == "gencode" {
        core::Gencode::new(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str(),
            None,
            None,
            None,
            None
        )
    } else {
        panic!("Unsupported annotation source: {}", reference_gene_annotation_source);
    };
    let scoring_method: caller::ReferenceTranscriptScoringMethod = reference_transcript_scoring_method.as_str().parse().unwrap();
    let selection_strategy: caller::ReferenceTranscriptSelectionStrategy = reference_transcript_selection_strategy.as_str().parse().unwrap();
    let tms: caller::TranscriptModelSet = caller::identify_variant_transcripts(
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
        num_threads,
        chunk_size,
        temp_dir.as_str()
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_assembled_transcripts: DataFrame = caller::assembled_transcript_records_to_dataframe(
                caller::build_assembled_transcript_records(&tms)
            );
            let df_exons: DataFrame = caller::exon_records_to_dataframe(
                caller::build_exon_records(&tms)
            );
            let df_introns: DataFrame = caller::intron_records_to_dataframe(
                caller::build_intron_records(&tms)
            );
            let df_reference_transcript_matches: DataFrame = caller::reference_transcript_match_records_to_dataframe(
                caller::build_reference_transcript_match_records(&tms)
            );
            let df_read_filter_status: DataFrame = caller::read_filter_status_records_to_dataframe(
                caller::build_read_filter_status_records(&tms)
            );
            let df_transcript_model_structures: DataFrame = caller::transcript_model_structure_records_to_dataframe(
                caller::build_transcript_model_structure_records(&tms)
            );
            let df_rna_variants: DataFrame = caller::rna_variant_records_to_dataframe(
                caller::build_rna_variant_records(&tms)
            );

            Ok((PyDataFrame(df_assembled_transcripts),
                PyDataFrame(df_exons),
                PyDataFrame(df_introns),
                PyDataFrame(df_reference_transcript_matches),
                PyDataFrame(df_read_filter_status),
                PyDataFrame(df_transcript_model_structures),
                PyDataFrame(df_rna_variants)))
        }
        "file" => {
            let assembled_transcripts_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcripts.tsv", output_prefix));
            let exons_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_exons.tsv", output_prefix));
            let introns_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_introns.tsv", output_prefix));
            let reference_transcript_matches_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_reference_transcript_matches.tsv", output_prefix));
            let read_filter_status_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_read_filter_status.tsv", output_prefix));
            let transcript_model_structures_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_transcript_model_structures.tsv", output_prefix));
            let rna_variants_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_rna_variants.tsv", output_prefix));

            core::write_tsv_file(
                caller::build_assembled_transcript_records(&tms),
                &assembled_transcripts_tsv_file
            );

            core::write_tsv_file(
                caller::build_exon_records(&tms),
                &exons_tsv_file
            );

            core::write_tsv_file(
                caller::build_intron_records(&tms),
                &introns_tsv_file
            );

            core::write_tsv_file(
                caller::build_reference_transcript_match_records(&tms),
                &reference_transcript_matches_tsv_file
            );

            core::write_tsv_file(
                caller::build_read_filter_status_records(&tms),
                &read_filter_status_tsv_file
            );

            core::write_tsv_file(
                caller::build_transcript_model_structure_records(&tms),
                &transcript_model_structures_tsv_file
            );

            core::write_tsv_file(
                caller::build_rna_variant_records(&tms),
                &rna_variants_tsv_file
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
