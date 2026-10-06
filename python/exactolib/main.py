# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement Exacto's main APIs.
"""


import gc
import pandas as pd
import polars as pl
from exactolib import exactolibrs
from typing import List, Tuple

from .constants import *
from .default import *
from .index.reference import build_reference_kmer_index
from .logging import get_logger
from .utilities import str2bool
from .variant_calling.peptide import identify_peptide_variants as identify_peptide_variants_


logger = get_logger(__name__)


def annotate_variant_calls(
        tsv_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        gene_types: List[str],
        gene_levels: List[int],
        transcript_types: List[str],
        transcript_levels: List[int],
        output_tsv_file: str,
        num_threads: int = ANNOTATE_VARS_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> pd.DataFrame:
    df_variants = exactolibrs.annotate_variant_calls(
        tsv_file=tsv_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        gene_types=gene_types,
        gene_levels=gene_levels,
        transcript_types=transcript_types,
        transcript_levels=transcript_levels,
        output_tsv_file=output_tsv_file,
        num_threads=num_threads,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def build_genome_variation_graph(
        df_variants: pl.DataFrame,
        fasta_file: str,
        output_fasta_file: str,
        sequence_prefix: str,
        remove_unknown_bases: bool,
        only_variant_sequences: bool,
        graph_type: str,
        num_threads: int,
        output_type: OutputType = OutputType.FILE,
        verbose: bool = True
):
    if output_type == OutputType.DATAFRAME:
        df_sequences = exactolibrs.build_genome_variation_graph(
            df_variants=df_variants,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            sequence_prefix=sequence_prefix,
            remove_unknown_bases=remove_unknown_bases,
            only_variant_sequences=only_variant_sequences,
            graph_type=graph_type,
            num_threads=num_threads,
            output_type=str(output_type),
            verbose=verbose
        )
        return df_sequences.to_pandas()
    if output_type == OutputType.VECTOR:
        sequences = exactolibrs.build_genome_variation_graph(
            df_variants=df_variants,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            sequence_prefix=sequence_prefix,
            remove_unknown_bases=remove_unknown_bases,
            only_variant_sequences=only_variant_sequences,
            graph_type=graph_type,
            num_threads=num_threads,
            output_type=str(output_type),
            verbose=verbose
        )
        return sequences
    if output_type == OutputType.FILE:
        exactolibrs.build_genome_variation_graph(
            df_variants=df_variants,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            sequence_prefix=sequence_prefix,
            remove_unknown_bases=remove_unknown_bases,
            only_variant_sequences=only_variant_sequences,
            graph_type=graph_type,
            num_threads=num_threads,
            output_type=str(output_type),
            verbose=verbose
        )
        return None


def build_transcriptome_variation_graph(
        df_transcript_structures: pl.DataFrame,
        fasta_file: str,
        output_fasta_file: str,
        graph_type: str,
        num_threads: int,
        output_type: OutputType = OutputType.FILE,
        batch_size: int = BUILD_TRANSCRIPTOME_VAR_GRAPH_BATCH_SIZE,
        verbose: bool = True
):
    if output_type == OutputType.DATAFRAME:
        df_sequences = exactolibrs.build_transcriptome_variation_graph(
            df_transcript_structures=df_transcript_structures,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            graph_type=graph_type,
            num_threads=num_threads,
            batch_size=batch_size,
            output_type=str(output_type),
            verbose=verbose
        )
        return df_sequences.to_pandas()
    if output_type == OutputType.VECTOR:
        sequences = exactolibrs.build_transcriptome_variation_graph(
            df_transcript_structures=df_transcript_structures,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            graph_type=graph_type,
            num_threads=num_threads,
            batch_size=batch_size,
            output_type=str(output_type),
            verbose=verbose
        )
        return sequences
    if output_type == OutputType.FILE:
        exactolibrs.build_transcriptome_variation_graph(
            df_transcript_structures=df_transcript_structures,
            fasta_file=fasta_file,
            output_fasta_file=output_fasta_file,
            graph_type=graph_type,
            num_threads=num_threads,
            batch_size=batch_size,
            output_type=str(output_type),
            verbose=verbose
        )
        return None


def cluster_rna_reads(
        bam_file: str,
        bai_file: str,
        reference_genome_fasta_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        output_dir: str,
        output_prefix: str,
        analyte_type: str = CLUSTER_RNA_READS_ANALYTE_TYPE,
        remove_unspliced_rnas: bool = str2bool(CLUSTER_RNA_READS_REMOVE_UNSPLICED_RNAS),
        max_batch_reads: int = CLUSTER_RNA_READS_MAX_BATCH_READS,
        max_locus_gap: int = CLUSTER_RNA_READS_MAX_LOCUS_GAP,
        unspliced_bin_size: int = CLUSTER_RNA_READS_UNSPLICED_BIN_SIZE,
        min_mapping_quality: int = CLUSTER_RNA_READS_MIN_MAPPING_QUALITY,
        min_terminal_softclip_length: int = CLUSTER_RNA_READS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        min_reads_per_cluster: int = CLUSTER_RNA_READS_MIN_READS_PER_CLUSTER,
        max_records: int = CLUSTER_RNA_READS_MAX_RECORDS,
        expected_sequencing_error: float = CLUSTER_RNA_READS_EXPECTED_SEQUENCING_ERROR,
        expected_slippage_probability: float = CLUSTER_RNA_READS_EXPECTED_SLIPPAGE_PROBABILITY,
        max_fpr: float = CLUSTER_RNA_READS_MAX_FPR,
        max_slippage_repeat_length: int = CLUSTER_RNA_READS_MAX_SLIPPAGE_REPEAT_LENGTH,
        min_size_proportion: float = CLUSTER_RNA_READS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CLUSTER_RNA_READS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance: int = CLUSTER_RNA_READS_MAX_INTRACHROMOSOMAL_DISTANCE,
        dna_variants_tsv_file: str = "",
        allowed_variants_tsv_file: str = "",
        allowed_variant_max_distance: int = CLUSTER_RNA_READS_ALLOWED_VARIANT_MAX_DISTANCE,
        mec_max_k: int = CLUSTER_RNA_READS_MEC_MAX_K,
        mec_num_restarts: int = CLUSTER_RNA_READS_MEC_NUM_RESTARTS,
        mec_max_iter: int = CLUSTER_RNA_READS_MEC_MAX_ITER,
        mec_seed: int = CLUSTER_RNA_READS_MEC_SEED,
        ts_min_homology: int = CLUSTER_RNA_READS_TS_MIN_HOMOLOGY,
        ts_soft_min_homology: int = CLUSTER_RNA_READS_TS_SOFT_MIN_HOMOLOGY,
        ts_max_breakpoint_dispersion: int = CLUSTER_RNA_READS_TS_MAX_BREAKPOINT_DISPERSION,
        max_unspliced_locus_gap: int = CLUSTER_RNA_READS_MAX_UNSPLICED_LOCUS_GAP,
        chunk_size: int = CLUSTER_RNA_READS_CHUNK_SIZE,
        soft_clip_removal: bool = str2bool(CLUSTER_RNA_READS_SOFT_CLIP_REMOVAL),
        soft_clip_max_boundary_distance: int = CLUSTER_RNA_READS_SOFT_CLIP_MAX_BOUNDARY_DISTANCE,
        soft_clip_bases_per_edit: int = CLUSTER_RNA_READS_SOFT_CLIP_BASES_PER_EDIT,
        soft_clip_min_partner_bases: int = CLUSTER_RNA_READS_SOFT_CLIP_MIN_PARTNER_BASES,
        unplaced_tail_removal: bool = str2bool(CLUSTER_RNA_READS_UNPLACED_TAIL_REMOVAL),
        unplaced_tail_min_ins_len: int = CLUSTER_RNA_READS_UNPLACED_TAIL_MIN_INS_LEN,
        bkpt_rescue: bool = str2bool(CLUSTER_RNA_READS_BKPT_RESCUE),
        bkpt_rescue_min_ins_len: int = CLUSTER_RNA_READS_BKPT_RESCUE_MIN_INS_LEN,
        bkpt_rescue_realignment_gap_open_score: int = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        bkpt_rescue_realignment_gap_extend_score: int = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        bkpt_rescue_realignment_k: int = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_K,
        bkpt_rescue_realignment_band_width: int = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        bkpt_rescue_realignment_min_score_fraction: float = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        bkpt_rescue_realignment_min_query_coverage: float = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        bkpt_rescue_realignment_min_placed_fraction: float = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION,
        bkpt_rescue_realignment_max_pieces: int = CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES,
        poa_match_score: int = CLUSTER_RNA_READS_POA_MATCH_SCORE,
        poa_mismatch_score: int = CLUSTER_RNA_READS_POA_MISMATCH_SCORE,
        poa_gap_open_score: int = CLUSTER_RNA_READS_POA_GAP_OPEN_SCORE,
        poa_gap_extend_score: int = CLUSTER_RNA_READS_POA_GAP_EXTEND_SCORE,
        min_reads: int = CLUSTER_RNA_READS_MIN_READS,
        min_total_depth: int = CLUSTER_RNA_READS_MIN_TOTAL_DEPTH,
        min_homopolymer_len: int = CLUSTER_RNA_READS_MIN_HOMOPOLYMER_LEN,
        min_dinucleotide_context_len: int = CLUSTER_RNA_READS_MIN_DINUCLEOTIDE_CONTEXT_LEN,
        template_switch_flank: int = CLUSTER_RNA_READS_TEMPLATE_SWITCH_FLANK,
        template_switch_foldback_min_stem: int = CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MIN_STEM,
        template_switch_foldback_max_loop_len: int = CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_LOOP_LEN,
        template_switch_foldback_max_distance: int = CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_DISTANCE,
        template_switch_foldback_slack: int = CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_SLACK,
        num_threads: int = CLUSTER_RNA_READS_NUM_THREADS,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    (df_clusters, df_summary, df_ref_transcripts, df_splice_junctions, df_variants,
     df_failed_variants, df_template_switch) = exactolibrs.cluster_rna_reads(
        bam_file=bam_file,
        bai_file=bai_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        output_dir=output_dir,
        output_prefix=output_prefix,
        analyte_type=analyte_type,
        remove_unspliced_rnas=remove_unspliced_rnas,
        max_batch_reads=max_batch_reads,
        max_locus_gap=max_locus_gap,
        unspliced_bin_size=unspliced_bin_size,
        min_mapping_quality=min_mapping_quality,
        min_terminal_softclip_length=min_terminal_softclip_length,
        min_reads_per_cluster=min_reads_per_cluster,
        max_records=max_records,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_fpr=max_fpr,
        max_slippage_repeat_length=max_slippage_repeat_length,
        min_size_proportion=min_size_proportion,
        max_ins_norm_edit_distance=max_ins_norm_edit_distance,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        dna_variants_tsv_file=dna_variants_tsv_file,
        allowed_variants_tsv_file=allowed_variants_tsv_file,
        allowed_variant_max_distance=allowed_variant_max_distance,
        mec_max_k=mec_max_k,
        mec_num_restarts=mec_num_restarts,
        mec_max_iter=mec_max_iter,
        mec_seed=mec_seed,
        ts_min_homology=ts_min_homology,
        ts_soft_min_homology=ts_soft_min_homology,
        ts_max_breakpoint_dispersion=ts_max_breakpoint_dispersion,
        max_unspliced_locus_gap=max_unspliced_locus_gap,
        chunk_size=chunk_size,
        soft_clip_removal=soft_clip_removal,
        soft_clip_max_boundary_distance=soft_clip_max_boundary_distance,
        soft_clip_bases_per_edit=soft_clip_bases_per_edit,
        soft_clip_min_partner_bases=soft_clip_min_partner_bases,
        unplaced_tail_removal=unplaced_tail_removal,
        unplaced_tail_min_ins_len=unplaced_tail_min_ins_len,
        bkpt_rescue=bkpt_rescue,
        bkpt_rescue_min_ins_len=bkpt_rescue_min_ins_len,
        bkpt_rescue_realignment_gap_open_score=bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_placed_fraction=bkpt_rescue_realignment_min_placed_fraction,
        bkpt_rescue_realignment_max_pieces=bkpt_rescue_realignment_max_pieces,
        poa_match_score=poa_match_score,
        poa_mismatch_score=poa_mismatch_score,
        poa_gap_open_score=poa_gap_open_score,
        poa_gap_extend_score=poa_gap_extend_score,
        min_reads=min_reads,
        min_total_depth=min_total_depth,
        min_homopolymer_len=min_homopolymer_len,
        min_dinucleotide_context_len=min_dinucleotide_context_len,
        template_switch_flank=template_switch_flank,
        template_switch_foldback_min_stem=template_switch_foldback_min_stem,
        template_switch_foldback_max_loop_len=template_switch_foldback_max_loop_len,
        template_switch_foldback_max_distance=template_switch_foldback_max_distance,
        template_switch_foldback_slack=template_switch_foldback_slack,
        num_threads=num_threads,
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_clusters, df_summary, df_ref_transcripts, df_splice_junctions, df_variants, df_failed_variants, df_template_switch


def correct_rna_reads(
        bam_file: str,
        clusters_tsv_file: str,
        cluster_reference_transcripts_tsv_file: str,
        cluster_splice_junctions_tsv_file: str,
        cluster_variants_tsv_file: str,
        output_dir: str,
        output_prefix: str,
        reference_gene_annotation_file: str = "",
        reference_gene_annotation_source: GeneAnnotationSource = GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly: str = "",
        reference_gene_annotation_version: str = "",
        max_correctable_event_len: int = CORRECT_RNA_READS_MAX_CORRECTABLE_EVENT_LEN,
        corrected_base_quality: int = CORRECT_RNA_READS_CORRECTED_BASE_QUALITY,
        trim_transcript_ends: bool = CORRECT_RNA_READS_TRIM_TRANSCRIPT_ENDS,
        num_threads: int = CORRECT_RNA_READS_NUM_THREADS,
        chunk_size: int = CORRECT_RNA_READS_CHUNK_SIZE
) -> str:
    """
    Error-correct RNA reads against the cluster each was assigned to.

    Args:
        bam_file                                :   BAM file (the alignment cluster-rna-reads ran on).
        clusters_tsv_file                       :   cluster-rna-reads clusters TSV file.
        cluster_reference_transcripts_tsv_file  :   cluster-rna-reads reference transcripts TSV file.
        cluster_splice_junctions_tsv_file       :   cluster-rna-reads splice junctions TSV file.
        cluster_variants_tsv_file               :   cluster-rna-reads passed variants TSV file.
        output_dir                              :   Output directory.
        output_prefix                           :   Output prefix.
        reference_gene_annotation_file          :   Reference gene annotation file ("" = none).
                                                    Read only with trim_transcript_ends.
        reference_gene_annotation_source        :   Reference gene annotation source.
        reference_gene_annotation_assembly      :   Reference gene annotation assembly.
        reference_gene_annotation_version       :   Reference gene annotation version.
        max_correctable_event_len               :   Uncalled insertions and deletions longer than this
                                                    are left as sequenced.
        corrected_base_quality                  :   Base quality assigned to a rewritten base.
        trim_transcript_ends                    :   Cut a read end that runs past an internal
                                                    annotated exon back to the exon.
        num_threads                             :   Number of threads.
        chunk_size                              :   Number of reads corrected per chunk.

    Returns:
        Path of the corrected reads FASTQ file (BGZF),
        {output_dir}/{output_prefix}_exacto_rna_corrected_reads.fastq.gz.
    """
    corrected_reads_fastq_file = exactolibrs.correct_rna_reads(
        bam_file=bam_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        clusters_tsv_file=clusters_tsv_file,
        cluster_reference_transcripts_tsv_file=cluster_reference_transcripts_tsv_file,
        cluster_splice_junctions_tsv_file=cluster_splice_junctions_tsv_file,
        cluster_variants_tsv_file=cluster_variants_tsv_file,
        output_dir=output_dir,
        output_prefix=output_prefix,
        max_correctable_event_len=max_correctable_event_len,
        corrected_base_quality=corrected_base_quality,
        trim_transcript_ends=trim_transcript_ends,
        num_threads=num_threads,
        chunk_size=chunk_size
    )
    return corrected_reads_fastq_file


def determine_rna_consensus(
        tsv_file: str,
        fastq_file: str,
        output_dir: str,
        output_prefix: str,
        match_score: int = DETERMINE_RNA_CONSENSUS_SCORE_MATCH,
        mismatch_score: int = DETERMINE_RNA_CONSENSUS_SCORE_MISMATCH,
        gap_open_score: int = DETERMINE_RNA_CONSENSUS_SCORE_GAP_OPEN,
        gap_extend_score: int = DETERMINE_RNA_CONSENSUS_SCORE_GAP_EXTEND,
        max_reads_per_cluster: int = DETERMINE_RNA_CONSENSUS_MAX_READS_PER_CLUSTER,
        seed: int = DETERMINE_RNA_CONSENSUS_SEED,
        orientation_kmer_size: int = DETERMINE_RNA_CONSENSUS_ORIENTATION_KMER_SIZE,
        num_threads: int = DETERMINE_RNA_CONSENSUS_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> pd.DataFrame:
    df_consensus = exactolibrs.determine_rna_consensus(
        tsv_file=tsv_file,
        fastq_file=fastq_file,
        output_dir=output_dir,
        output_prefix=output_prefix,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
        max_reads_per_cluster=max_reads_per_cluster,
        seed=seed,
        orientation_kmer_size=orientation_kmer_size,
        num_threads=num_threads,
        output_type=str(output_type)
    )
    return df_consensus


def identify_somatic_dna_variants(
        bam_file: str,
        bam_bai_file: str,
        control_bam_files: List[str],
        control_bam_bai_files: List[str],
        fasta_file: str,
        output_tsv_file: str,
        regions: List[Tuple[str, int, int]],
        min_reads: int = CALL_SOMATIC_DNA_VARS_MIN_READS,
        min_mapping_quality: int = CALL_SOMATIC_DNA_VARS_MIN_MAPPING_QUALITY,
        min_base_quality: int = CALL_SOMATIC_DNA_VARS_MIN_BASE_QUALITY,
        min_total_depth: int = CALL_SOMATIC_DNA_VARS_MIN_TOTAL_DEPTH,
        min_alt_allele_fraction: float = CALL_SOMATIC_DNA_VARS_MIN_ALT_ALLELE_FRACTION,
        min_size_proportion: float = CALL_SOMATIC_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_SOMATIC_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance: int = CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        max_slippage_repeat_length: int = CALL_SOMATIC_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH,
        min_terminal_softclip_length: int = CALL_SOMATIC_DNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        bkpt_rescue: bool = str2bool(CALL_SOMATIC_DNA_VARS_BKPT_RESCUE),
        bkpt_rescue_min_ins_len: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MIN_INS_LEN,
        bkpt_rescue_max_ins_len: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MAX_INS_LEN,
        bkpt_rescue_search_distance: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_SEARCH_DISTANCE,
        bkpt_rescue_realignment_gap_open_score: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        bkpt_rescue_realignment_gap_extend_score: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        bkpt_rescue_realignment_k: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_K,
        bkpt_rescue_realignment_band_width: int = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        bkpt_rescue_realignment_min_score_fraction: float = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        bkpt_rescue_realignment_min_query_coverage: float = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        bkpt_rescue_realignment_min_span_proportion: float = CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SPAN_PROPORTION,
        poa_match_score: int = CALL_SOMATIC_DNA_VARS_POA_MATCH_SCORE,
        poa_mismatch_score: int = CALL_SOMATIC_DNA_VARS_POA_MISMATCH_SCORE,
        poa_gap_open_score: int = CALL_SOMATIC_DNA_VARS_POA_GAP_OPEN_SCORE,
        poa_gap_extend_score: int = CALL_SOMATIC_DNA_VARS_POA_GAP_EXTEND_SCORE,
        read_depth_max_merge_distance: int = CALL_SOMATIC_DNA_VARS_READ_DEPTH_MAX_MERGE_DISTANCE,
        min_homopolymer_len: int = CALL_SOMATIC_DNA_VARS_MIN_HOMOPOLYMER_LEN,
        min_dinucleotide_context_len: int = CALL_SOMATIC_DNA_VARS_MIN_DINUCLEOTIDE_CONTEXT_LEN,
        max_control_reads: int = CALL_SOMATIC_DNA_VARS_MAX_CONTROL_READS,
        num_threads: int = CALL_SOMATIC_DNA_VARS_NUM_THREADS,
        chunk_size: int = CALL_SOMATIC_DNA_VARS_CHUNK_SIZE,
        max_records: int = CALL_SOMATIC_DNA_VARS_MAX_RECORDS,
        expected_variant_allele_fraction: float = CALL_SOMATIC_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION,
        expected_mutation_rate: float = CALL_SOMATIC_DNA_VARS_EXPECTED_MUTATION_RATE,
        expected_sequencing_error: float = CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
        expected_slippage_probability: float = CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
        max_fpr: float = CALL_SOMATIC_DNA_VARS_MAX_FPR,
        temp_dir: str = "",
        apply_infinite_sites_assumption: bool = str2bool(CALL_SOMATIC_DNA_VARS_INFINITE_SITES_ASSUMPTION),
        output_type: OutputType = OutputType.FILE
):
    """
    Identify case-specific DNA variants.

    Args:
        bam_file                            :   Case BAM file.
        bam_bai_file                        :   Case BAM.BAI file.
        control_bam_files                   :   List of control BAM files.
        control_bam_bai_files               :   List of control BAM.BAI files.
        output_tsv_file                     :   Output TSV file.
        chromosomes                         :   A list of chromosomes to scan for the presence of DNA variants (default: None).
                                                If left unspecified (i.e. None), all chromosomes in the BAM file will be considered.
        min_reads                           :   Minimum number of reads.
        min_mapping_quality                 :   Minimum mapping quality.
        min_average_base_quality            :   Minimum average base quality.
        min_size_proportion                 :   Minimum size proportion.
        max_ins_norm_edit_distance          :   Maximum insertion edit distance.
        max_intrachromosomal_distance       :   Maximum intrachromosomal distance for clustering variants.
        apply_infinite_sites_assumption     :   If True, any variant in the case BAM file that shares the same position
                                                as a variant in any of the control BAM file will be filtered out.
        num_threads                         :   Number of threads.
        temp_dir                            :   Temp directory (default: TMPDIR).
        output_type                         :   Output type ('file' or 'dataframe').

    Returns:
        If output_type is 'dataframe', then a Pandas DataFrame.
    """
    assert len(control_bam_files) > 0, 'control_bam_files cannot be empty.'
    assert len(control_bam_files) == len(control_bam_bai_files), 'len(control_bam_files) must be equal to len(control_bam_bai_files).'
    df_variants = exactolibrs.identify_somatic_dna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        control_bam_files=control_bam_files,
        control_bam_bai_files=control_bam_bai_files,
        fasta_file=fasta_file,
        output_tsv_file=output_tsv_file,
        regions=regions,
        min_reads=min_reads,
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_total_depth=min_total_depth,
        min_alt_allele_fraction=min_alt_allele_fraction,
        min_size_proportion=min_size_proportion,
        max_ins_norm_edit_distance=max_ins_norm_edit_distance,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_slippage_repeat_length=max_slippage_repeat_length,
        min_terminal_softclip_length=min_terminal_softclip_length,
        bkpt_rescue=bkpt_rescue,
        bkpt_rescue_min_ins_len=bkpt_rescue_min_ins_len,
        bkpt_rescue_max_ins_len=bkpt_rescue_max_ins_len,
        bkpt_rescue_search_distance=bkpt_rescue_search_distance,
        bkpt_rescue_realignment_gap_open_score=bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_span_proportion=bkpt_rescue_realignment_min_span_proportion,
        poa_match_score=poa_match_score,
        poa_mismatch_score=poa_mismatch_score,
        poa_gap_open_score=poa_gap_open_score,
        poa_gap_extend_score=poa_gap_extend_score,
        read_depth_max_merge_distance=read_depth_max_merge_distance,
        min_homopolymer_len=min_homopolymer_len,
        min_dinucleotide_context_len=min_dinucleotide_context_len,
        max_control_reads=max_control_reads,
        num_threads=num_threads,
        chunk_size=chunk_size,
        max_records=max_records,
        expected_variant_allele_fraction=expected_variant_allele_fraction,
        expected_mutation_rate=expected_mutation_rate,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_fpr=max_fpr,
        apply_infinite_sites_assumption=apply_infinite_sites_assumption,
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def identify_germline_dna_variants(
        bam_file: str,
        bam_bai_file: str,
        fasta_file: str,
        output_tsv_file: str,
        regions: List[Tuple[str,int,int]],
        min_reads: int = CALL_GERMLINE_DNA_VARS_MIN_READS,
        min_mapping_quality: int = CALL_GERMLINE_DNA_VARS_MIN_MAPPING_QUALITY,
        min_base_quality: int = CALL_GERMLINE_DNA_VARS_MIN_BASE_QUALITY,
        min_total_depth: int = CALL_GERMLINE_DNA_VARS_MIN_TOTAL_DEPTH,
        min_alt_allele_fraction: float = CALL_GERMLINE_DNA_VARS_MIN_ALT_ALLELE_FRACTION,
        min_size_proportion: float = CALL_GERMLINE_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_GERMLINE_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance: int = CALL_GERMLINE_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        max_slippage_repeat_length: int = CALL_GERMLINE_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH,
        min_terminal_softclip_length: int = CALL_GERMLINE_DNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        bkpt_rescue: bool = str2bool(CALL_GERMLINE_DNA_VARS_BKPT_RESCUE),
        bkpt_rescue_min_ins_len: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_MIN_INS_LEN,
        bkpt_rescue_max_ins_len: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_MAX_INS_LEN,
        bkpt_rescue_search_distance: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_SEARCH_DISTANCE,
        bkpt_rescue_realignment_gap_open_score: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        bkpt_rescue_realignment_gap_extend_score: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        bkpt_rescue_realignment_k: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_K,
        bkpt_rescue_realignment_band_width: int = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        bkpt_rescue_realignment_min_score_fraction: float = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        bkpt_rescue_realignment_min_query_coverage: float = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        bkpt_rescue_realignment_min_span_proportion: float = CALL_GERMLINE_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SPAN_PROPORTION,
        poa_match_score: int = CALL_GERMLINE_DNA_VARS_POA_MATCH_SCORE,
        poa_mismatch_score: int = CALL_GERMLINE_DNA_VARS_POA_MISMATCH_SCORE,
        poa_gap_open_score: int = CALL_GERMLINE_DNA_VARS_POA_GAP_OPEN_SCORE,
        poa_gap_extend_score: int = CALL_GERMLINE_DNA_VARS_POA_GAP_EXTEND_SCORE,
        read_depth_max_merge_distance: int = CALL_GERMLINE_DNA_VARS_READ_DEPTH_MAX_MERGE_DISTANCE,
        min_homopolymer_len: int = CALL_GERMLINE_DNA_VARS_MIN_HOMOPOLYMER_LEN,
        min_dinucleotide_context_len: int = CALL_GERMLINE_DNA_VARS_MIN_DINUCLEOTIDE_CONTEXT_LEN,
        num_threads: int = CALL_GERMLINE_DNA_VARS_NUM_THREADS,
        chunk_size: int = CALL_GERMLINE_DNA_VARS_CHUNK_SIZE,
        max_records: int = CALL_GERMLINE_DNA_VARS_MAX_RECORDS,
        expected_variant_allele_fraction: float = CALL_GERMLINE_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION,
        expected_mutation_rate: float = CALL_GERMLINE_DNA_VARS_EXPECTED_MUTATION_RATE,
        expected_sequencing_error: float = CALL_GERMLINE_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
        expected_slippage_probability: float = CALL_GERMLINE_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
        max_fpr: float = CALL_GERMLINE_DNA_VARS_MAX_FPR,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
) -> pd.DataFrame:
    """
    Identify DNA variants.

    Args:
        bam_file                            :   BAM file.
        bam_bai_file                        :   BAM.BAI file.
        output_tsv_file                     :   Output TSV file.
        chromosomes                         :   A list of chromosomes to scan for the presence of DNA variants (default: None).
                                                If left unspecified (i.e. None), all chromosomes in the BAM file will be considered.
        min_reads                           :   Minimum number of supporting reads.
        min_mapping_quality                 :   Minimum mapping quality.
        min_average_base_quality            :   Minimum average base quality.
        min_size_proportion                 :   Minimum size proportion.
        max_ins_norm_edit_distance          :   Maximum insertion edit distance.
        max_intrachromosomal_distance       :   Maximum intrachromosomal distance for clustering variants.
        num_threads                         :   Number of threads.
        chromosomes                         :   Chromosomes in which to identify variants (default: []).
                                                If left unspecified, all chromosomes in the BAM file will be considered.
        temp_dir                            :   Temp directory (default: TMPDIR).
        output_type                         :   Output type ('file' or 'dataframe').

    Returns:
        If output_type is 'dataframe', then a Pandas DataFrame.
    """
    df_variants = exactolibrs.identify_germline_dna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        fasta_file=fasta_file,
        output_tsv_file=output_tsv_file,
        regions=regions,
        min_reads=min_reads,
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_total_depth=min_total_depth,
        min_alt_allele_fraction=min_alt_allele_fraction,
        min_size_proportion=min_size_proportion,
        max_ins_norm_edit_distance=max_ins_norm_edit_distance,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_slippage_repeat_length=max_slippage_repeat_length,
        min_terminal_softclip_length=min_terminal_softclip_length,
        bkpt_rescue=bkpt_rescue,
        bkpt_rescue_min_ins_len=bkpt_rescue_min_ins_len,
        bkpt_rescue_max_ins_len=bkpt_rescue_max_ins_len,
        bkpt_rescue_search_distance=bkpt_rescue_search_distance,
        bkpt_rescue_realignment_gap_open_score=bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_span_proportion=bkpt_rescue_realignment_min_span_proportion,
        poa_match_score=poa_match_score,
        poa_mismatch_score=poa_mismatch_score,
        poa_gap_open_score=poa_gap_open_score,
        poa_gap_extend_score=poa_gap_extend_score,
        read_depth_max_merge_distance=read_depth_max_merge_distance,
        min_homopolymer_len=min_homopolymer_len,
        min_dinucleotide_context_len=min_dinucleotide_context_len,
        num_threads=num_threads,
        chunk_size=chunk_size,
        max_records=max_records,
        expected_variant_allele_fraction=expected_variant_allele_fraction,
        expected_mutation_rate=expected_mutation_rate,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_fpr=max_fpr,
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def identify_rna_transcript_variants(
        bam_file: str,
        reference_genome_fasta_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        output_dir: str,
        output_prefix: str,
        # cluster_variants_tsv_file: str = "",
        chunk_size: int = CALL_RNA_VARS_CHUNK_SIZE,
        min_mapping_quality: int = CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        min_terminal_softclip_length: int = CALL_RNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        bkpt_rescue: bool = str2bool(CALL_RNA_VARS_BKPT_RESCUE),
        bkpt_rescue_min_ins_len: int = CALL_RNA_VARS_BKPT_RESCUE_MIN_INS_LEN,
        bkpt_rescue_realignment_gap_open_score: int = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        bkpt_rescue_realignment_gap_extend_score: int = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        bkpt_rescue_realignment_k: int = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_K,
        bkpt_rescue_realignment_band_width: int = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        bkpt_rescue_realignment_min_score_fraction: float = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        bkpt_rescue_realignment_min_query_coverage: float = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        bkpt_rescue_realignment_min_placed_fraction: float = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION,
        bkpt_rescue_realignment_max_pieces: int = CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES,
        predict_nonsense_mediated_decay: bool = str2bool(CALL_RNA_VARS_PREDICT_NONSENSE_MEDIATED_DECAY),
        translation_strategy: str = CALL_RNA_VARS_TRANSLATION_STRATEGY,
        start_codons: List[str] = CALL_RNA_VARS_START_CODONS,
        nmd_distance_threshold: int = CALL_RNA_VARS_NMD_DISTANCE_THRESHOLD,
        dna_variants_tsv_files: List[str] = [],
        dna_variant_match_buffer: int = CALL_RNA_VARS_DNA_VARIANT_MATCH_BUFFER,
        num_threads: int = CALL_RNA_VARS_NUM_THREADS,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Identify RNA transcript variants.

    Args:
        bam_file                            :   BAM file of assembled transcripts.
        reference_genome_fasta_file         :   Reference genome FASTA file.
        reference_gene_annotation_file      :   Reference gene annotation file.
        reference_gene_annotation_source    :   Reference gene annotation source.
        reference_gene_annotation_assembly  :   Reference gene annotation assembly.
        reference_gene_annotation_version   :   Reference gene annotation version.
        output_dir                          :   Output directory.
        output_prefix                       :   Output prefix.
        cluster_variants_tsv_file           :   Passed cluster variants TSV file written by
                                                'cluster_rna_reads'. When given, a transcript keeps
                                                only the variants its own cluster called, so every
                                                read name in the BAM file must be a cluster ID
                                                (default: "" keeps every variant).
        chunk_size                          :   Chunk size for variant calling.
        min_mapping_quality                 :   Minimum mapping quality.
        min_terminal_softclip_length        :   Terminal soft-clip runs shorter than this emit no variant record.
        dna_variants_tsv_files              :   DNA variant TSV files written by 'call_germline_dna_vars'
                                                or 'call_somatic_dna_vars' (column 'origin' required).
                                                An RNA variant with the same chromosomes, operations,
                                                variant type and sequence as a DNA variant, each position
                                                within 'dna_variant_match_buffer' bases (strand not
                                                compared), carries the DNA variant's origin in the RNA
                                                variants' 'origin' column (default: [] leaves it empty).
        dna_variant_match_buffer            :   Bases each position of a matching DNA variant may differ by.
        num_threads                         :   Number of threads.
        temp_dir                            :   Temp directory (default: TMPDIR).
        output_type                         :   Output type ('file' or 'dataframe').

    Returns:
        If output_type is 'dataframe', then
        Pandas DataFrame of assembled transcripts,
        Pandas DataFrame of exons,
        Pandas DataFrame of introns,
        Pandas DataFrame of reference transcript matches,
        Pandas DataFrame of read filter status,
        Pandas DataFrame of transcript model structures,
        Pandas DataFrame of RNA variants,
        Pandas DataFrame of nonsense-mediated decay predictions (one row per open reading frame)
    """
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants,
     df_nmd_predictions) = exactolibrs.identify_rna_transcript_variants(
        bam_file=bam_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        # cluster_variants_tsv_file=cluster_variants_tsv_file,
        output_dir=output_dir,
        output_prefix=output_prefix,
        min_mapping_quality=min_mapping_quality,
        min_terminal_softclip_length=min_terminal_softclip_length,
        bkpt_rescue=bkpt_rescue,
        bkpt_rescue_min_ins_len=bkpt_rescue_min_ins_len,
        bkpt_rescue_realignment_gap_open_score=bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_placed_fraction=bkpt_rescue_realignment_min_placed_fraction,
        bkpt_rescue_realignment_max_pieces=bkpt_rescue_realignment_max_pieces,
        predict_nonsense_mediated_decay=predict_nonsense_mediated_decay,
        translation_strategy=str(translation_strategy),
        start_codons=start_codons,
        nmd_distance_threshold=nmd_distance_threshold,
        dna_variants_tsv_files=dna_variants_tsv_files,
        dna_variant_match_buffer=dna_variant_match_buffer,
        num_threads=num_threads,
        temp_dir=temp_dir,
        output_type=str(output_type),
        chunk_size=chunk_size
    )
    return (df_assembled_transcripts.to_pandas(),
            df_exons.to_pandas(),
            df_introns.to_pandas(),
            df_reference_transcript_matches.to_pandas(),
            df_read_filter_status.to_pandas(),
            df_transcript_model_structures.to_pandas(),
            df_rna_variants.to_pandas(),
            df_nmd_predictions.to_pandas())


def identify_peptide_variants(
        proteoforms_tsv_file: str,
        reference_proteome_fasta_file: str,
        min_k: int,
        max_k: int,
        num_processes: int
) -> pd.DataFrame:
    """
    Identify mutant peptide k-mers from a proteoforms TSV.

    Reads the `ProteoformRecord` TSV produced by `translate_transcripts`,
    extracts every k-mer (min_k..max_k) that overlaps at least one mutant
    amino-acid index, and filters out k-mers that occur in the reference
    proteome.

    Args:
        proteoforms_tsv_file        :   Path to the proteoforms TSV written by
                                        translate-transcripts. Must contain at
                                        minimum:
                                            proteoform_id,
                                            amino_acid_sequence,
                                            mutant_amino_acid_intervals,
                                            assembled_transcript_variant_ids,
                                            dna_variant_ids.
        reference_proteome_fasta_file: Path to the reference proteome FASTA.
        min_k                       :   Minimum peptide length.
        max_k                       :   Maximum peptide length.
        num_processes               :   Number of worker processes.

    Returns:
        pd.DataFrame with columns:
            mutant_peptide_id, proteoform_id, mutant_peptide_sequence, k,
            amino_acid_index_start, amino_acid_index_end,
            assembled_transcript_variant_ids, dna_variant_ids.
    """
    reference_kmer_set = build_reference_kmer_index(
        reference_fasta_file=reference_proteome_fasta_file,
        min_k=min_k,
        max_k=max_k,
        num_processes=num_processes
    )
    # keep_default_na=False preserves empty strings as "" instead of NaN,
    # so the worker doesn't have to special-case NaN for the variant-ID columns.
    df_proteoforms = pd.read_csv(
        proteoforms_tsv_file,
        sep='\t',
        low_memory=False,
        memory_map=True,
        keep_default_na=False
    )
    df_mutant_peptides = identify_peptide_variants_(
        df_proteoforms=df_proteoforms,
        reference_kmer_set=reference_kmer_set,
        min_k=min_k,
        max_k=max_k,
        num_processes=num_processes
    )
    return df_mutant_peptides


def integrate_variants(
        dna_variants_tsv_file: str,
        rna_variants_tsv_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        output_tsv_file: str,
        max_exon_offset: int = INTEGRATE_VARS_MAX_EXON_OFFSET,
        max_transcript_boundary_offset: int = INTEGRATE_VARS_MAX_TRANSCRIPT_BOUNDARY_OFFSET,
        max_intergenic_distance: int = INTEGRATE_VARS_MAX_INTERGENIC_DISTANCE,
        num_threads: int = INTEGRATE_VARS_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> pd.DataFrame:
    """
    Integrate DNA and RNA variants.

    Args:
        annotated_dna_variant_callset_tsv_file      :   Annotated DNA variant callset TSV file.
        rna_variant_callset_tsv_file                :   RNA variant callset TSV file.
        reference_gene_annotation_file              :   Reference gene annotation TSV file.
        reference_gene_annotation_source            :   Reference gene annotation TSV file source.
        output_tsv_file                             :   Output TSV file.
        max_exon_offset                             :   Maximum exon offset.
        max_transcript_boundary_offset              :   Maximum transcript boundary offset.
        max_intergenic_distance                     :   Maximum intergenic distance.
        num_threads                                 :   Number of threads.
        output_type                                 :   Output type ('file' or 'dataframe').

    Returns:
        Pandas DataFrame with the following columns:
            'rna_variant_call_id,
            'dna_variant_call_id',
            'distance',
            'rna_variant_position',
            'dna_variant_position'
    """
    df_integration = exactolibrs.integrate_dna_rna_variants(
        dna_variants_tsv_file=dna_variants_tsv_file,
        rna_variants_tsv_file=rna_variants_tsv_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        output_tsv_file=output_tsv_file,
        max_exon_offset=max_exon_offset,
        max_transcript_boundary_offset=max_transcript_boundary_offset,
        max_intergenic_distance=max_intergenic_distance,
        num_threads=num_threads,
        output_type=str(output_type)
    )
    return df_integration.to_pandas()


def quantify_rna_abundances(
        clusters_tsv_file: str,
        cluster_reference_matches_tsv_file: str,
        cluster_splice_junctions_tsv_file: str,
        cluster_variants_tsv_file: str,
        consensus_rna_reference_matches_tsv_file: str,
        output_dir: str,
        output_prefix: str,
        pseudo_count: float = QUANTIFY_RNA_ABUNDANCES_PSEUDO_COUNT,
        max_iter: int = QUANTIFY_RNA_ABUNDANCES_MAX_ITER,
        tol: float = QUANTIFY_RNA_ABUNDANCES_TOL,
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df_abundances, df_abundances_ref_txs = exactolibrs.quantify_rna_abundances(
        clusters_tsv_file=clusters_tsv_file,
        cluster_reference_matches_tsv_file=cluster_reference_matches_tsv_file,
        cluster_splice_junctions_tsv_file=cluster_splice_junctions_tsv_file,
        cluster_variants_tsv_file=cluster_variants_tsv_file,
        consensus_rna_reference_matches_tsv_file=consensus_rna_reference_matches_tsv_file,
        output_dir=output_dir,
        output_prefix=output_prefix,
        pseudo_count=pseudo_count,
        max_iter=max_iter,
        tol=tol,
        output_type=str(output_type)
    )
    return (df_abundances.to_pandas(),
            df_abundances_ref_txs.to_pandas())


def remove_unspliced_rnas(
        bam_file: str,
        bam_bai_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        gene_types: List[str],
        gene_levels: List[int],
        transcript_types: List[str],
        transcript_levels: List[int],
        output_bam_file: str,
        output_bam_bai_file: str,
        num_threads: int,
        min_mapping_quality: int,
        write_output_bam_file: bool = True
) -> List[str]:
    read_names_to_keep = exactolibrs.remove_unspliced_rnas(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        gene_types=gene_types,
        gene_levels=gene_levels,
        transcript_types=transcript_types,
        transcript_levels=transcript_levels,
        output_bam_file=output_bam_file,
        output_bam_bai_file=output_bam_bai_file,
        num_threads=num_threads,
        min_mapping_quality=min_mapping_quality,
        write_output_bam_file=write_output_bam_file
    )
    return read_names_to_keep


def stitch_reference_transcripts(
        bam_file: str,
        reference_genome_fasta_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        output_dir: str,
        output_prefix: str,
        min_mapping_quality: int = STITCH_REFERENCE_TRANSCRIPTS_MIN_MAPPING_QUALITY,
        min_num_splice_junction_matches: int = STITCH_REFERENCE_TRANSCRIPTS_MIN_NUM_SPLICE_JUNCTION_MATCHES,
        pas_search_size: int = STITCH_REFERENCE_TRANSCRIPTS_PAS_SEARCH_SIZE,
        pas_hexamers: List[str] = STITCH_REFERENCE_TRANSCRIPTS_PAS_HEXAMERS,
        pas_start_offset_range: Tuple[int, int] = STITCH_REFERENCE_TRANSCRIPTS_PAS_START_OFFSET_RANGE,
        polya_tail_min_adenosine_fraction: float = STITCH_REFERENCE_TRANSCRIPTS_POLYA_TAIL_MIN_ADENOSINE_FRACTION,
        polya_window_size: int = STITCH_REFERENCE_TRANSCRIPTS_POLYA_WINDOW_SIZE,
        polya_min_consecutive_adenosine: int = STITCH_REFERENCE_TRANSCRIPTS_POLYA_MIN_CONSECUTIVE_ADENOSINE,
        num_threads: int = STITCH_REFERENCE_TRANSCRIPTS_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> pd.DataFrame:
    """
    Stitch reference transcript sequence onto the degraded ends of aligned RNA sequences.

    Every mapped read in the BAM file is one transcript, identified by its read name. In
    'file' mode the stitched transcripts TSV and FASTA are written; the stitched sequences no
    longer match the read coordinates of any call-rna-transcript-vars output computed on the
    same BAM file: realign the stitched FASTA and re-run call-rna-transcript-vars before
    translate-transcripts.

    Args:
        bam_file                            :   BAM file of the aligned RNA sequences (e.g. realigned consensus sequences).
        reference_genome_fasta_file         :   Reference genome FASTA file.
        reference_gene_annotation_file      :   Reference gene annotation file.
        reference_gene_annotation_source    :   Reference gene annotation source.
        reference_gene_annotation_assembly  :   Reference gene annotation assembly.
        reference_gene_annotation_version   :   Reference gene annotation version.
        output_dir                          :   Output directory.
        output_prefix                       :   Output prefix.
        min_mapping_quality                 :   Reads mapped below this are left out.
        min_num_splice_junction_matches     :   Matches below this are treated as no match.
        pas_search_size                     :   Read bases searched for a polyadenylation signal hexamer.
        pas_hexamers                        :   Polyadenylation signal hexamers.
        pas_start_offset_range              :   Inclusive (min, max) distance of a hexamer's start from the searched end.
        polya_tail_min_adenosine_fraction   :   Minimum A fraction of the read bases past the 3' terminal aligned base.
        polya_window_size                   :   Genomic bases scanned on each side of an intronic 3' end for an A-tract.
        polya_min_consecutive_adenosine     :   A-run length in that window that calls internal priming.
        num_threads                         :   Number of threads.
        output_type                         :   Output type ('file' or 'dataframe').

    Returns:
        Pandas DataFrame of stitched transcripts, one row per read. Empty when output_type is 'file'.
    """
    df_stitched = exactolibrs.stitch_reference_transcripts(
        bam_file=bam_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        min_mapping_quality=min_mapping_quality,
        min_num_splice_junction_matches=min_num_splice_junction_matches,
        pas_search_size=pas_search_size,
        pas_hexamers=pas_hexamers,
        pas_start_offset_range=tuple(pas_start_offset_range),
        polya_tail_min_adenosine_fraction=polya_tail_min_adenosine_fraction,
        polya_window_size=polya_window_size,
        polya_min_consecutive_adenosine=polya_min_consecutive_adenosine,
        num_threads=num_threads,
        output_dir=output_dir,
        output_prefix=output_prefix,
        output_type=str(output_type)
    )
    return df_stitched.to_pandas()


def translate_fastx_file(
        fastx_file: str,
        output_fasta_file: str,
        output_tsv_file: str,
        strategy: TranslationStrategy = TranslationStrategy(TRANSLATE_STRATEGY),
        start_codons: List[str] = ["AUG"],
        num_threads: int = TRANSLATE_NUM_THREADS
):
    """
    Translate a long-read RNA-seq FASTX file into peptide sequences.

    Args:
        fastx_file          :   FASTX file.
        output_fasta_file   :   Output FASTX file.
        output_tsv_file     :   Output TSV file.
        strategy            :   Translation strategy.
        start_codons        :   List of start codons (default: ["AUG"]).
        num_threads         :   Number of threads.
    """
    tsv_file = exactolibrs.translate_fastx_file(
        fastx_file=fastx_file,
        output_fasta_file=output_fasta_file,
        output_tsv_file=output_tsv_file,
        strategy=str(strategy),
        start_codons=start_codons,
        num_threads=num_threads
    )
    gc.collect()


def translate_sequence(
        rna_sequence: str,
        strategy: TRANSLATE_STRATEGY,
        start_codons: List[str] = ["AUG"]
) -> List[Tuple[str, int, int]]:
    """
    Translate a single RNA sequence.

    Args:
        rna_sequence    :   RNA sequence.
        strategy        :   Translation strategy.
        start_codons    :   List of start codons (default: ["AUG"]).

    Returns:
        List of (peptide_sequence, orf_start, orf_end) tuples — one per
        primary structure produced by the strategy.
    """
    translations = exactolibrs.translate_sequence(
        rna_sequence=rna_sequence,
        strategy=str(strategy),
        start_codons=start_codons
    )
    return translations


def translate_transcripts(
        assembled_transcript_model_alignments_tsv_file: str,
        assembled_transcript_variants_tsv_file: str,
        dna_variants_tsv_file: str,
        integrated_variants_tsv_file: str,
        strategy: TRANSLATE_STRATEGY,
        output_dir: str,
        output_prefix: str,
        assembled_transcript_support_tsv_file: str = "",
        rna_consensus_tsv_file: str = "",
        stitched_transcripts_tsv_file: str = "",
        start_codons: List[str] = TRANSLATE_START_CODONS,
        num_threads: int = TRANSLATE_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Translate transcript structures.

    Transcript sequences and their supporting read names come from exactly one of
    assembled_transcript_support_tsv_file or rna_consensus_tsv_file.

    Args:
        assembled_transcript_support_tsv_file           :   Assembled transcript support TSV file
                                                            (externally assembled transcripts).
        rna_consensus_tsv_file                          :   RNA consensus TSV file (output from
                                                            exacto determine-rna-consensus).
        stitched_transcripts_tsv_file                   :   Stitched transcripts TSV file written by
                                                            exacto stitch-reference-transcripts (optional). Its stitched
                                                            sequences replace the transcript sequences, so the
                                                            call-rna-transcript-vars TSV files must come from the
                                                            realigned stitched FASTA. Adds reference-stitched
                                                            provenance columns to both outputs.
        assembled_transcript_model_alignments_tsv_file  :   Assembled transcript model alignments TSV file.
        assembled_transcript_variants_tsv_file          :   Assembled transcript variants TSV file.
        dna_variants_tsv_file                           :   DNA variant records TSV file.
        integrated_variants_tsv_file                    :   Integrated variant records TSV file.
        strategy                                        :   Translation strategy.
        output_tsv_file                                 :   Output TSV file path (written in 'file' mode).
        output_fasta_file                               :   Output FASTA file path (empty = skip; only honored in 'file' mode).
        start_codons                                    :   List of start codons (default: ["AUG"]).
        num_threads                                     :   Number of threads.
        output_type                                     :   'file' writes TSV+FASTA and returns an empty DataFrame;
                                                    '       dataframe' returns the populated DataFrame and writes no files.

    Returns:
        Pandas DataFrame. Empty in 'file' mode; populated with one row per
        primary structure in 'dataframe' mode.
    """
    if bool(assembled_transcript_support_tsv_file) == bool(rna_consensus_tsv_file):
        raise ValueError(
            "Pass exactly one of assembled_transcript_support_tsv_file (externally assembled "
            "transcripts) or rna_consensus_tsv_file (output of determine-rna-consensus)."
        )

    df_proteoforms, df_proteoform_nucleotides = exactolibrs.translate_transcripts(
        assembled_transcript_support_tsv_file=assembled_transcript_support_tsv_file,
        rna_consensus_tsv_file=rna_consensus_tsv_file,
        stitched_transcripts_tsv_file=stitched_transcripts_tsv_file,
        assembled_transcript_model_alignments_tsv_file=assembled_transcript_model_alignments_tsv_file,
        assembled_transcript_variants_tsv_file=assembled_transcript_variants_tsv_file,
        dna_variants_tsv_file=dna_variants_tsv_file,
        integrated_variants_tsv_file=integrated_variants_tsv_file,
        strategy=str(strategy),
        start_codons=start_codons,
        output_dir=output_dir,
        output_prefix=output_prefix,
        num_threads=num_threads,
        output_type=str(output_type)
    )

    return df_proteoforms.to_pandas(), df_proteoform_nucleotides.to_pandas()

