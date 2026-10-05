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
use exacto::cluster::prelude as cluster;
use exacto::core::prelude as core;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use crate::functions::{io_error, resolve_temp_dir};


/// `cluster_variants_tsv_file` is optional: the passed-variants table `cluster-rna-reads`
/// wrote for the clusters the BAM's transcripts were built from. When given, a transcript
/// keeps only the variants its own cluster called, so every read name in the BAM must be a
/// cluster id. An empty string turns the filter off.
#[pyfunction]
pub fn identify_rna_transcript_variants(
    py: Python,
    bam_file: String,
    reference_genome_fasta_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
//     cluster_variants_tsv_file: String,
    output_dir: String,
    output_prefix: String,
    min_mapping_quality: u16,
    min_terminal_softclip_length: u32,
    bkpt_rescue: bool,
    bkpt_rescue_min_ins_len: u32,
    bkpt_rescue_realignment_gap_open_score: i32,
    bkpt_rescue_realignment_gap_extend_score: i32,
    bkpt_rescue_realignment_k: u32,
    bkpt_rescue_realignment_band_width: u32,
    bkpt_rescue_realignment_min_score_fraction: f64,
    bkpt_rescue_realignment_min_query_coverage: f64,
    bkpt_rescue_realignment_min_placed_fraction: f64,
    bkpt_rescue_realignment_max_pieces: u32,
    predict_nonsense_mediated_decay: bool,
    translation_strategy: String,
    start_codons: Vec<String>,
    nmd_distance_threshold: u32,
    dna_variants_tsv_files: Vec<String>,
    dna_variant_match_buffer: u32,
    num_threads: usize,
    temp_dir: String,
    output_type: String,
    chunk_size: usize
) -> PyResult<(PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame)> {
    let _ = py;
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
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported annotation source: {}", reference_gene_annotation_source
        )));
    };

    let Ok(translation_strategy) = core::TranslationStrategy::from_str(translation_strategy.as_str()) else {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported translation strategy: {:?}. Expected \"longest_orf\" or \"all_orfs\".",
            translation_strategy
        )));
    };

    // Every option field is a parameter, so the Python defaults are the only defaults.
    let options: caller::IdentifyRNATranscriptVariantsOptions = caller::IdentifyRNATranscriptVariantsOptions {
        calling: caller::RNAVariantCallingOptions {
            chunk_size: chunk_size,
            min_mapping_quality: min_mapping_quality,
            min_terminal_soft_clip_ins_len: min_terminal_softclip_length,
            bkpt_rescue: bkpt_rescue,
            bkpt_rescue_min_ins_len: bkpt_rescue_min_ins_len,
            bkpt_rescue_realignment_gap_open_score: bkpt_rescue_realignment_gap_open_score,
            bkpt_rescue_realignment_gap_extend_score: bkpt_rescue_realignment_gap_extend_score,
            bkpt_rescue_realignment_k: bkpt_rescue_realignment_k,
            bkpt_rescue_realignment_band_width: bkpt_rescue_realignment_band_width,
            bkpt_rescue_realignment_min_score_fraction: bkpt_rescue_realignment_min_score_fraction,
            bkpt_rescue_realignment_min_query_coverage: bkpt_rescue_realignment_min_query_coverage,
            bkpt_rescue_realignment_min_placed_fraction: bkpt_rescue_realignment_min_placed_fraction,
            bkpt_rescue_realignment_max_pieces: bkpt_rescue_realignment_max_pieces
        },
        annotation: caller::RNAAnnotationOptions {
            predict_nonsense_mediated_decay: predict_nonsense_mediated_decay,
            translation_strategy: translation_strategy,
            start_codons: start_codons.into_iter().collect(),
            nmd_distance_threshold: nmd_distance_threshold
        }
    };

//     // The cluster table spells chromosomes by name; the pipeline compares graph operations by
//     // the BAM's chromosome ids, so the operations are rebuilt on those ids here.
//     let cluster_rna_variants: Option<HashMap<usize, HashSet<caller::GraphOperation>>> = if cluster_variants_tsv_file.is_empty() {
//         None
//     } else {
//         let chromosome_names_map = core::create_chromosome_names_map(bam_file.as_str());
//         let mut operations_by_cluster: HashMap<usize, HashSet<caller::GraphOperation>> = HashMap::new();
//         for record in cluster::load_rna_read_cluster_variant_records(cluster_variants_tsv_file.as_str()) {
//             let mut chromosome_ids: Vec<u16> = Vec::with_capacity(2);
//             for chromosome_name in [&record.chromosome_1, &record.chromosome_2] {
//                 let Some(chromosome_id) = chromosome_names_map.get_by_left(chromosome_name) else {
//                     return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
//                         "Chromosome {} in {} is not in the header of {}",
//                         chromosome_name, cluster_variants_tsv_file, bam_file
//                     )));
//                 };
//                 chromosome_ids.push(*chromosome_id);
//             }
//             let parse_error = |column: &str, value: &str| -> PyErr {
//                 PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
//                     "Unsupported {} in {}: {}", column, cluster_variants_tsv_file, value
//                 ))
//             };
//             let operation: caller::GraphOperation = caller::GraphOperation::new(
//                 chromosome_ids[0],
//                 record.position_1,
//                 core::Strand::from_str(&record.strand_1)
//                     .map_err(|_| parse_error("strand_1", &record.strand_1))?,
//                 caller::GraphOperationType::from_str(&record.operation_type_1)
//                     .map_err(|_| parse_error("operation_type_1", &record.operation_type_1))?,
//                 chromosome_ids[1],
//                 record.position_2,
//                 core::Strand::from_str(&record.strand_2)
//                     .map_err(|_| parse_error("strand_2", &record.strand_2))?,
//                 caller::GraphOperationType::from_str(&record.operation_type_2)
//                     .map_err(|_| parse_error("operation_type_2", &record.operation_type_2))?,
//                 record.sequence.clone(),
//                 caller::VariantType::from_str(&record.variant_type)
//                     .map_err(|_| parse_error("variant_type", &record.variant_type))?
//             );
//             operations_by_cluster
//                 .entry(record.cluster_id)
//                 .or_default()
//                 .insert(operation);
//         }
//         Some(operations_by_cluster)
//     };

    // Optional DNA variants: output TSVs of `call-germline-dna-vars` / `call-somatic-dna-vars`.
    // An RNA variant they match carries their origins in the variants table.
    let dna_variant_records: Vec<caller::DNAVariantRecord> = dna_variants_tsv_files
        .iter()
        .flat_map(|dna_variants_tsv_file| caller::load_dna_variant_records(dna_variants_tsv_file.as_str()))
        .collect();

    let tms: caller::TranscriptModelSet = caller::identify_rna_transcript_variants(
        bam_file.as_str(),
        reference_genome_fasta_file.as_str(),
        &gene_annotator,
//         cluster_rna_variants.as_ref(),
        &options,
        num_threads,
        resolve_temp_dir(temp_dir.as_str()).as_str()
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_assembled_transcripts: DataFrame = caller::assembled_transcript_records_to_dataframe(
                caller::build_assembled_transcript_records(&tms)
            );
            let df_assembled_transcript_exons: DataFrame = caller::assembled_transcript_exon_records_to_dataframe(
                caller::build_assembled_transcript_exon_records(&tms)
            );
            let df_assembled_transcript_splice_junctions: DataFrame = caller::assembled_transcript_splice_junction_records_to_dataframe(
                caller::build_assembled_transcript_splice_junction_records(&tms)
            );
            let df_assembled_transcript_reference_transcript_matches: DataFrame = caller::assembled_transcript_reference_transcript_match_records_to_dataframe(
                caller::build_assembled_transcript_reference_transcript_match_records(&tms)
            );
            let df_assembled_transcript_filter_status: DataFrame = caller::assembled_transcript_filter_status_records_to_dataframe(
                caller::build_assembled_transcript_filter_status_records(&tms)
            );
            let df_assembled_transcript_nmd_predictions: DataFrame = caller::assembled_transcript_nonsense_mediated_decay_records_to_dataframe(
                caller::build_assembled_transcript_nonsense_mediated_decay_records(&tms)
            );
            let df_assembled_transcript_model_alignments: DataFrame = caller::assembled_transcript_model_alignment_records_to_dataframe(
                caller::build_assembled_transcript_alignment_records(&tms)
            );
            let df_assembled_transcript_variants: DataFrame = caller::assembled_transcript_variant_records_to_dataframe(
                caller::annotate_dna_variant_origins(
                    caller::build_assembled_transcript_variant_records(
                        &tms,
                        Some(reference_genome_fasta_file.as_str()),
                        20
                    ),
                    &dna_variant_records,
                    dna_variant_match_buffer
                )
            );

            Ok((PyDataFrame(df_assembled_transcripts),
                PyDataFrame(df_assembled_transcript_exons),
                PyDataFrame(df_assembled_transcript_splice_junctions),
                PyDataFrame(df_assembled_transcript_reference_transcript_matches),
                PyDataFrame(df_assembled_transcript_filter_status),
                PyDataFrame(df_assembled_transcript_model_alignments),
                PyDataFrame(df_assembled_transcript_variants),
                PyDataFrame(df_assembled_transcript_nmd_predictions)))
        }
        "file" => {
            let assembled_transcripts_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcripts.tsv", output_prefix));
            let assembled_transcript_exons_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_exons.tsv", output_prefix));
            let assembled_transcript_splice_junctions_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_splice_junctions.tsv", output_prefix));
            let assembled_transcript_reference_transcript_matches_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_reference_transcript_matches.tsv", output_prefix));
            let assembled_transcript_filter_status_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_filter_status.tsv", output_prefix));
            let assembled_transcript_model_alignments_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_model_alignments.tsv", output_prefix));
            let assembled_transcript_variants_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_variants.tsv", output_prefix));
            let assembled_transcript_nmd_predictions_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_assembled_transcript_nmd_predictions.tsv", output_prefix));

            core::write_tsv_table(
                caller::build_assembled_transcript_records(&tms),
                &assembled_transcripts_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_exon_records(&tms),
                &assembled_transcript_exons_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_splice_junction_records(&tms),
                &assembled_transcript_splice_junctions_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_reference_transcript_match_records(&tms),
                &assembled_transcript_reference_transcript_matches_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_filter_status_records(&tms),
                &assembled_transcript_filter_status_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_alignment_records(&tms),
                &assembled_transcript_model_alignments_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::annotate_dna_variant_origins(
                    caller::build_assembled_transcript_variant_records(
                        &tms,
                        Some(reference_genome_fasta_file.as_str()),
                        20
                    ),
                    &dna_variant_records,
                    dna_variant_match_buffer
                ),
                &assembled_transcript_variants_tsv_file
            ).map_err(io_error)?;

            core::write_tsv_table(
                caller::build_assembled_transcript_nonsense_mediated_decay_records(&tms),
                &assembled_transcript_nmd_predictions_tsv_file
            ).map_err(io_error)?;

            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
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
