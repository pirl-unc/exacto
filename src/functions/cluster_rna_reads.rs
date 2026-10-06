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
use crate::functions::resolve_temp_dir;


#[pyfunction]
pub fn cluster_rna_reads(
    py: Python,
    bam_file: String,
    bai_file: String,
    reference_genome_fasta_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    output_dir: String,
    output_prefix: String,
    analyte_type: String,
    remove_unspliced_rnas: bool,
    max_batch_reads: usize,
    max_locus_gap: u32,
    unspliced_bin_size: u32,
    min_mapping_quality: u16,
    min_terminal_softclip_length: u32,
    min_reads_per_cluster: usize,
    expected_sequencing_error: f64,
    expected_slippage_probability: f64,
    max_records: usize,
    max_fpr: f64,
    max_slippage_repeat_length: u32,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_intrachromosomal_distance: u32,
    dna_variants_tsv_file: String,
    allowed_variants_tsv_file: String,
    allowed_variant_max_distance: u32,
    mec_max_k: usize,
    mec_num_restarts: usize,
    mec_max_iter: usize,
    mec_seed: usize,
    ts_min_homology: u32,
    ts_soft_min_homology: u32,
    ts_max_breakpoint_dispersion: u32,
    max_unspliced_locus_gap: u32,
    chunk_size: usize,
    soft_clip_removal: bool,
    soft_clip_max_boundary_distance: u32,
    soft_clip_bases_per_edit: u32,
    soft_clip_min_partner_bases: u32,
    unplaced_tail_removal: bool,
    unplaced_tail_min_ins_len: u32,
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
    poa_match_score: i32,
    poa_mismatch_score: i32,
    poa_gap_open_score: i32,
    poa_gap_extend_score: i32,
    min_reads: usize,
    min_total_depth: usize,
    min_homopolymer_len: u32,
    min_dinucleotide_context_len: u32,
    template_switch_flank: u32,
    template_switch_foldback_min_stem: u32,
    template_switch_foldback_max_loop_len: u32,
    template_switch_foldback_max_distance: u32,
    template_switch_foldback_slack: u32,
    num_threads: usize,
    temp_dir: String,
    output_type: String
) -> PyResult<(PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame, PyDataFrame)> {
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

    // Optional DNA variant keep list: an output TSV of `call-germline-dna-vars` or
    // `call-somatic-dna-vars`. Empty string means no list.
    let dna_variant_records: Option<Vec<caller::DNAVariantRecord>> =
        if dna_variants_tsv_file.as_str() == "" {
            None
        } else {
            Some(caller::load_dna_variant_records(dna_variants_tsv_file.as_str()))
        };

    // Optional allow list: the `_variants_passed.tsv` of a pass-1 `cluster-rna-reads` run, which
    // restricts a pass-2 run to the variants pass 1 called. Empty string means no list.
    let allowed_variant_records: Option<Vec<cluster::RNAReadClusterVariantRecord>> =
        if allowed_variants_tsv_file.as_str() == "" {
            None
        } else {
            Some(cluster::load_rna_read_cluster_variant_records(allowed_variants_tsv_file.as_str()))
        };

    let analyte_type: core::AnalyteType = match analyte_type.as_str() {
        "cdna" => core::AnalyteType::CDNA,
        "rna" => core::AnalyteType::RNA,
        other => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unsupported value for analyte_type: {} (cdna or rna)", other
            )));
        }
    };

    let options: cluster::ClusterRNAReadsOptions = cluster::ClusterRNAReadsOptions {
        analyte_type: analyte_type,
        seed: mec_seed,
        remove_unspliced_rnas: remove_unspliced_rnas,
        max_batch_reads: max_batch_reads,
        max_locus_gap: max_locus_gap,
        unspliced_bin_size: unspliced_bin_size,
        error_model: caller::SequencingErrorModel {
            sequencing_error: expected_sequencing_error,
            slippage_prob: expected_slippage_probability
        },
        junction: cluster::JunctionClusteringOptions {
            min_reads: min_reads_per_cluster,
            max_unspliced_locus_gap: max_unspliced_locus_gap
        },
        calling: cluster::RNAVariantCallingOptions {
            chunk_size: chunk_size,
            max_records: max_records,
            min_mapping_quality: min_mapping_quality,
            max_clustering_distance: max_intrachromosomal_distance,
            max_ins_norm_edit_distance: max_ins_norm_edit_distance,
            min_size_proportion: min_size_proportion,
            min_terminal_soft_clip_ins_len: min_terminal_softclip_length,
            soft_clip_removal: soft_clip_removal,
            soft_clip_max_boundary_distance: soft_clip_max_boundary_distance,
            soft_clip_bases_per_edit: soft_clip_bases_per_edit,
            soft_clip_min_partner_bases: soft_clip_min_partner_bases,
            unplaced_tail_removal: unplaced_tail_removal,
            unplaced_tail_min_ins_len: unplaced_tail_min_ins_len,
            bkpt_rescue: bkpt_rescue,
            bkpt_rescue_min_ins_len: bkpt_rescue_min_ins_len,
            bkpt_rescue_realignment_gap_open_score: bkpt_rescue_realignment_gap_open_score,
            bkpt_rescue_realignment_gap_extend_score: bkpt_rescue_realignment_gap_extend_score,
            bkpt_rescue_realignment_k: bkpt_rescue_realignment_k,
            bkpt_rescue_realignment_band_width: bkpt_rescue_realignment_band_width,
            bkpt_rescue_realignment_min_score_fraction: bkpt_rescue_realignment_min_score_fraction,
            bkpt_rescue_realignment_min_query_coverage: bkpt_rescue_realignment_min_query_coverage,
            bkpt_rescue_realignment_min_placed_fraction: bkpt_rescue_realignment_min_placed_fraction,
            bkpt_rescue_realignment_max_pieces: bkpt_rescue_realignment_max_pieces,
            poa_match_score: poa_match_score,
            poa_mismatch_score: poa_mismatch_score,
            poa_gap_open_score: poa_gap_open_score,
            poa_gap_extend_score: poa_gap_extend_score
        },
        filtering: cluster::RNAVariantFilteringOptions {
            max_fpr: max_fpr,
            min_reads: min_reads,
            min_total_depth: min_total_depth,
            max_slippage_repeat_len: max_slippage_repeat_length,
            min_homopolymer_len: min_homopolymer_len,
            min_dinucleotide_context_len: min_dinucleotide_context_len,
            template_switch_flank: template_switch_flank,
            template_switch_junction_min_homology: ts_min_homology,
            template_switch_junction_soft_min_homology: ts_soft_min_homology,
            template_switch_junction_max_breakpoint_dispersion: ts_max_breakpoint_dispersion,
            template_switch_foldback_min_stem: template_switch_foldback_min_stem,
            template_switch_foldback_max_loop_len: template_switch_foldback_max_loop_len,
            template_switch_foldback_max_distance: template_switch_foldback_max_distance,
            template_switch_foldback_slack: template_switch_foldback_slack,
            allowed_variant_max_distance: allowed_variant_max_distance
        },
        phasing: cluster::PhasingOptions {
            mec_max_k: mec_max_k,
            mec_num_restart: mec_num_restarts,
            mec_max_iter: mec_max_iter
        }
    };

    let cluster_set: cluster::RNAReadClusterSet = cluster::cluster_rna_reads(
        bam_file.as_str(),
        bai_file.as_str(),
        reference_genome_fasta_file.as_str(),
        &gene_annotator,
        dna_variant_records.as_ref(),
        allowed_variant_records.as_ref(),
        &options,
        num_threads,
        resolve_temp_dir(temp_dir.as_str()).as_str(),
        output_dir.as_str(),
        output_prefix.as_str()
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_clusters: DataFrame = cluster::rna_read_cluster_id_records_to_dataframe(
                cluster::build_rna_read_cluster_id_records(&cluster_set)
            );
            let df_summary: DataFrame = cluster::rna_read_cluster_summary_records_to_dataframe(
                cluster::build_rna_read_cluster_summary_records(&cluster_set)
            );
            let df_reference_genes_transcripts: DataFrame = cluster::rna_read_cluster_reference_gene_transcript_records_to_dataframe(
                cluster::build_rna_read_cluster_reference_gene_transcript_records(&cluster_set)
            );
            let df_splice_junctions: DataFrame = cluster::rna_read_cluster_splice_junction_records_to_dataframe(
                cluster::build_rna_read_cluster_splice_junction_records(&cluster_set)
            );
            let df_variants: DataFrame = cluster::rna_read_cluster_variant_records_to_dataframe(
                cluster::build_rna_read_cluster_variant_records(&cluster_set)
            );
            let df_failed_variants: DataFrame = cluster::rna_read_cluster_failed_variant_records_to_dataframe(
                cluster::build_rna_read_cluster_failed_variant_records(&cluster_set)
            );
            let df_template_switch: DataFrame = cluster::rna_read_cluster_template_switch_records_to_dataframe(
                cluster::build_rna_read_cluster_template_switch_records(&cluster_set)
            );
            Ok((PyDataFrame(df_clusters),
                PyDataFrame(df_summary),
                PyDataFrame(df_reference_genes_transcripts),
                PyDataFrame(df_splice_junctions),
                PyDataFrame(df_variants),
                PyDataFrame(df_failed_variants),
                PyDataFrame(df_template_switch)))
        }
        "file" => {
            // `cluster_rna_reads` wrote the tables to `output_dir`.
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
