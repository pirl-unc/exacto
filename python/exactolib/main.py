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
import os
import pandas as pd
import polars as pl
from exactolib import exactolibrs
from typing import List, Tuple
from .constants import *
from .default import *
from .index.reference import build_reference_kmer_index
from .logging import get_logger
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
        min_base_quality: float = CALL_SOMATIC_DNA_VARS_MIN_BASE_QUALITY,
        min_total_depth: int = CALL_SOMATIC_DNA_VARS_MIN_TOTAL_DEPTH,
        min_alt_allele_fraction: float = CALL_SOMATIC_DNA_VARS_MIN_ALT_ALLELE_FRACTION,
        min_size_proportion: float = CALL_SOMATIC_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_SOMATIC_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance_tau: int = CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_intrachromosomal_distance: int = CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_interchromosomal_distance: int = CALL_SOMATIC_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        max_slippage_repeat_length: int = CALL_SOMATIC_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH,
        num_threads: int = CALL_SOMATIC_DNA_VARS_NUM_THREADS,
        chunk_size: int = CALL_SOMATIC_DNA_VARS_CHUNK_SIZE,
        max_records: int = CALL_SOMATIC_DNA_VARS_MAX_RECORDS,
        expected_variant_allele_fraction: float = CALL_SOMATIC_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION,
        expected_mutation_rate: float = CALL_SOMATIC_DNA_VARS_EXPECTED_MUTATION_RATE,
        expected_sequencing_error: float = CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
        expected_slippage_probability: float = CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
        max_f1_fraction: float = CALL_SOMATIC_DNA_VARS_MAX_F1_FRACTION,
        max_fpr: float = CALL_SOMATIC_DNA_VARS_MAX_FPR,
        temp_dir: str = "",
        apply_infinite_sites_assumption: bool = True,
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
        max_intrachromosomal_distance_tau   :   tau for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_intrachromosomal_distance       :   Maximum intrachromosomal distance.
                                                d_max for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_interchromosomal_distance       :   Maximum interchromosomal distance.
        apply_infinite_sites_assumption     :   If True, any variant in the case BAM file that shares the same position
                                                as a variant in any of the control BAM file will be filtered out.
        num_threads                         :   Number of threads.
        temp_dir                            :   Temp directory (default: TMPDIR).
        output_type                         :   Output type ('file' or 'dataframe').

    Returns:
        If output_type is 'dataframe', then a Pandas DataFrame.
    """
    assert(len(control_bam_files) > 0, 'control_bam_files cannot be empty.')
    assert(len(control_bam_files) == len(control_bam_bai_files), 'len(control_bam_files) must be equal to len(control_bam_bai_files).')
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
        max_intrachromosomal_distance_tau=max_intrachromosomal_distance_tau,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_interchromosomal_distance=max_interchromosomal_distance,
        max_slippage_repeat_length=max_slippage_repeat_length,
        num_threads=num_threads,
        chunk_size=chunk_size,
        max_records=max_records,
        expected_variant_allele_fraction=expected_variant_allele_fraction,
        expected_mutation_rate=expected_mutation_rate,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_f1_fraction=max_f1_fraction,
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
        min_base_quality: float = CALL_GERMLINE_DNA_VARS_MIN_BASE_QUALITY,
        min_total_depth: int = CALL_GERMLINE_DNA_VARS_MIN_TOTAL_DEPTH,
        min_alt_allele_fraction: float = CALL_GERMLINE_DNA_VARS_MIN_ALT_ALLELE_FRACTION,
        min_size_proportion: float = CALL_GERMLINE_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_GERMLINE_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance_tau: int = CALL_GERMLINE_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_intrachromosomal_distance: int = CALL_GERMLINE_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        max_interchromosomal_distance: int = CALL_GERMLINE_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        max_slippage_repeat_length: int = CALL_GERMLINE_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH,
        num_threads: int = CALL_GERMLINE_DNA_VARS_NUM_THREADS,
        chunk_size: int = CALL_GERMLINE_DNA_VARS_CHUNK_SIZE,
        max_records: int = CALL_GERMLINE_DNA_VARS_MAX_RECORDS,
        expected_variant_allele_fraction: float = CALL_GERMLINE_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION,
        expected_mutation_rate: float = CALL_GERMLINE_DNA_VARS_EXPECTED_MUTATION_RATE,
        expected_sequencing_error: float = CALL_GERMLINE_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
        expected_slippage_probability: float = CALL_GERMLINE_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
        max_f1_fraction: float = CALL_GERMLINE_DNA_VARS_MAX_F1_FRACTION,
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
        max_intrachromosomal_distance_tau   :   tau for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_intrachromosomal_distance       :   Maximum intrachromosomal distance.
                                                d_max for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_interchromosomal_distance       :   Maximum interchromosomal distance.
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
        max_intrachromosomal_distance_tau=max_intrachromosomal_distance_tau,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_interchromosomal_distance=max_interchromosomal_distance,
        max_slippage_repeat_length=max_slippage_repeat_length,
        num_threads=num_threads,
        chunk_size=chunk_size,
        max_records=max_records,
        expected_variant_allele_fraction=expected_variant_allele_fraction,
        expected_mutation_rate=expected_mutation_rate,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_f1_fraction=max_f1_fraction,
        max_fpr=max_fpr,
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def identify_rna_variants(
        bam_file: str,
        bam_bai_file: str,
        reference_genome_fasta_file: str,
        reference_gene_annotation_file: str,
        reference_gene_annotation_source: GeneAnnotationSource,
        reference_gene_annotation_assembly: str,
        reference_gene_annotation_version: str,
        output_dir: str,
        output_prefix: str,
        chunk_size: int = CALL_RNA_VARS_CHUNK_SIZE,
        reference_transcript_scoring_method: ReferenceTranscriptScoringMethod = ReferenceTranscriptScoringMethod(CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SCORING_METHOD),
        reference_transcript_selection_strategy: ReferenceTranscriptSelectionStrategy = ReferenceTranscriptSelectionStrategy(CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SELECTION_STRATEGY),
        reference_transcript_top_k: int = CALL_RNA_VARS_REFERENCE_TRANSCRIPT_TOP_K,
        reference_transcript_threshold: float = CALL_RNA_VARS_REFERENCE_TRANSCRIPT_THRESHOLD,
        min_mapping_quality: int = CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        min_average_base_quality: float = CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        num_threads: int = CALL_RNA_VARS_NUM_THREADS,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Identify RNA variants.

    Args:
        bam_file                            :   BAM file.
        bam_bai_file                        :   BAM.BAI file.
        genes_txt_file                      :   Genes TXT file.
        read_transcript_tsv_file            :   TSV file with the following columns:
                                                'read_id', 'transcript_id', 'probability'.
        min_reads                           :   Minimum number of reads.
        min_mapping_quality                 :   Minimum mapping quality.
        min_average_base_quality            :   Minimum average base quality.
        min_size_proportion                 :   Minimum size proportion.
        max_ins_norm_edit_distance          :   Maximum insertion edit distance.
        max_intrachromosomal_distance_tau   :   tau for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_intrachromosomal_distance       :   Maximum intrachromosomal distance.
                                                d_max for the clustering maximum distance formula:
                                                d = d_max * (1 - e^(-1*variant_size / tau))
        max_interchromosomal_distance       :   Maximum interchromosomal distance.
        num_threads                         :   Number of threads.
        chromosomes                         :   Chromosomes in which to identify variants (default: []).
                                                If left unspecified, all chromosomes in the BAM file will be considered.
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
        Pandas DataFrame of RNA variants
    """
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = exactolibrs.identify_rna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        output_dir=output_dir,
        output_prefix=output_prefix,
        reference_transcript_scoring_method=str(reference_transcript_scoring_method),
        reference_transcript_selection_strategy=str(reference_transcript_selection_strategy),
        reference_transcript_top_k=reference_transcript_top_k,
        reference_transcript_threshold=reference_transcript_threshold,
        min_mapping_quality=min_mapping_quality,
        min_average_base_quality=min_average_base_quality,
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
            df_rna_variants.to_pandas())


def identify_peptide_variants(
        primary_structures_tsv_file: str,
        reference_proteome_fasta_file: str,
        min_k: int,
        max_k: int,
        num_processes: int
) -> pd.DataFrame:
    """
    Identify mutant peptide k-mers from a primary-structures TSV.

    Reads the `PrimaryStructureRecord` TSV produced by `translate_structures`,
    extracts every k-mer (min_k..max_k) that overlaps at least one mutant
    amino-acid index, and filters out k-mers that occur in the reference
    proteome.

    Args:
        primary_structures_tsv_file :   Path to the primary-structures TSV.
                                        Must contain at minimum:
                                            primary_structure_id,
                                            amino_acid_sequence,
                                            mutant_amino_acid_intervals,
                                            rna_variant_ids,
                                            dna_variant_ids.
        reference_proteome_fasta_file: Path to the reference proteome FASTA.
        min_k                       :   Minimum peptide length.
        max_k                       :   Maximum peptide length.
        num_processes               :   Number of worker processes.

    Returns:
        pd.DataFrame with columns:
            mutant_peptide_id, primary_structure_id, mutant_peptide_sequence, k,
            amino_acid_index_start, amino_acid_index_end,
            rna_variant_ids, dna_variant_ids.
    """
    reference_kmer_set = build_reference_kmer_index(
        reference_fasta_file=reference_proteome_fasta_file,
        min_k=min_k,
        max_k=max_k,
        num_processes=num_processes
    )
    # keep_default_na=False preserves empty strings as "" instead of NaN,
    # so the worker doesn't have to special-case NaN for the variant-ID columns.
    df_primary_structures = pd.read_csv(
        primary_structures_tsv_file,
        sep='\t',
        low_memory=False,
        memory_map=True,
        keep_default_na=False
    )
    df_mutant_peptides = identify_peptide_variants_(
        df_primary_structures=df_primary_structures,
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
        num_threads: int = CALL_RNA_VARS_NUM_THREADS,
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


def remove_unspliced_rnas(
        bam_file: str,
        bam_bai_file: str,
        fasta_file: str,
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
        output_fasta_file: str,
        num_threads: int,
        min_mapping_quality: int
):
    exactolibrs.remove_unspliced_rnas(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        fasta_file=fasta_file,
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
        output_fasta_file=output_fasta_file,
        num_threads=num_threads,
        min_mapping_quality=min_mapping_quality
    )


def translate_fasta_file(
        fasta_file: str,
        temp_dir: str,
        strategy: TranslationStrategy = TranslationStrategy(TRANSLATE_STRATEGY),
        start_codons: List[str] = ["AUG"],
        num_threads: int = TRANSLATE_NUM_THREADS
) -> pd.DataFrame:
    """
    Translate a long-read RNA-seq FASTA file into peptide sequences.

    Args:
        fasta_file      :   FASTA file.
        strategy        :   Translation strategy.
        start_codons    :   List of start codons (default: ["AUG"]).
        num_threads     :   Number of threads.

    Returns:
        Pandas DataFrame with one row per `PrimaryStructure`, including
        columns: 'primary_structure_id', 'amino_acid_sequence',
        'transcript_id', 'assembled_transcript_name',
        'assembled_transcript_sequence', 'orf_start', 'orf_end',
        'reference_gene_names', 'reference_transcript_ids',
        'rna_variant_ids', 'rna_variants', 'dna_variant_ids',
        'dna_variants', 'rna_read_names', etc.
    """
    tsv_file = exactolibrs.translate_fasta_file(
        fasta_file=fasta_file,
        strategy=str(strategy),
        start_codons=start_codons,
        num_threads=num_threads,
        temp_dir=temp_dir
    )
    gc.collect()
    # Empty translation set → csv writer never emits the header row, so
    # the file is 0 bytes. Return an empty DataFrame in that case.
    if os.path.getsize(tsv_file) == 0:
        return pd.DataFrame()
    return pl.read_csv(tsv_file, separator='\t').to_pandas()


def translate_fastq_file(
        fastq_file: str,
        temp_dir: str,
        strategy: TranslationStrategy = TranslationStrategy(TRANSLATE_STRATEGY),
        start_codons: List[str] = ["AUG"],
        num_threads: int = TRANSLATE_NUM_THREADS
) -> pd.DataFrame:
    """
    Translate a long-read RNA-seq FASTQ file into peptide sequences.

    Args:
        fastq_file      :   FASTQ file (may be gzipped).
        strategy        :   Translation strategy.
        start_codons    :   List of start codons (default: ["AUG"]).
        num_threads     :   Number of threads.

    Returns:
        Pandas DataFrame with one row per `PrimaryStructure`, including
        columns: 'primary_structure_id', 'amino_acid_sequence',
        'transcript_id', 'assembled_transcript_name',
        'assembled_transcript_sequence', 'orf_start', 'orf_end',
        'reference_gene_names', 'reference_transcript_ids',
        'rna_variant_ids', 'rna_variants', 'dna_variant_ids',
        'dna_variants', 'rna_read_names', etc.
    """
    tsv_file = exactolibrs.translate_fastq_file(
        fastq_file=fastq_file,
        strategy=str(strategy),
        start_codons=start_codons,
        num_threads=num_threads,
        temp_dir=temp_dir
    )
    gc.collect()
    if os.path.getsize(tsv_file) == 0:
        return pd.DataFrame()
    return pl.read_csv(tsv_file, separator='\t').to_pandas()


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


def translate_structures(
        rna_assembly_support_tsv_file: str,
        transcript_model_structures_tsv_file: str,
        rna_variants_tsv_file: str,
        dna_variants_tsv_file: str,
        integrated_variants_tsv_file: str,
        strategy: TRANSLATE_STRATEGY,
        output_dir: str,
        output_prefix: str,
        start_codons: List[str] = ["AUG"],
        num_threads: int = TRANSLATE_NUM_THREADS,
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Translate transcript structures.

    Args:
        assembled_transcript_support_tsv_file   :   Assembled transcript support TSV file.
        transcript_model_structures_tsv_file    :   Transcript model structures TSV file.
        rna_variants_tsv_file                   :   RNA variant records TSV file.
        dna_variants_tsv_file                   :   DNA variant records TSV file.
        integrated_variants_tsv_file            :   Integrated variant records TSV file.
        strategy                                :   Translation strategy.
        output_tsv_file                         :   Output TSV file path (written in 'file' mode).
        output_fasta_file                       :   Output FASTA file path (empty = skip; only honored in 'file' mode).
        start_codons                            :   List of start codons (default: ["AUG"]).
        num_threads                             :   Number of threads.
        output_type                             :   'file' writes TSV+FASTA and returns an empty DataFrame;
                                                    'dataframe' returns the populated DataFrame and writes no files.

    Returns:
        Pandas DataFrame. Empty in 'file' mode; populated with one row per
        primary structure in 'dataframe' mode.
    """
    df_ps, df_ps_nucleotides = exactolibrs.translate_structures(
        rna_assembly_support_tsv_file=rna_assembly_support_tsv_file,
        transcript_model_structures_tsv_file=transcript_model_structures_tsv_file,
        rna_variants_tsv_file=rna_variants_tsv_file,
        dna_variants_tsv_file=dna_variants_tsv_file,
        integrated_variants_tsv_file=integrated_variants_tsv_file,
        strategy=str(strategy),
        start_codons=start_codons,
        output_dir=output_dir,
        output_prefix=output_prefix,
        num_threads=num_threads,
        output_type=str(output_type)
    )

    return df_ps.to_pandas(), df_ps_nucleotides.to_pandas()

