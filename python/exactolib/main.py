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
from typing import Dict, List, Optional, Tuple

from .constants import *
from .default import *
from .logging import get_logger
from .utilities import get_chromosomes,get_kmers


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
        temp_dir: str = "",
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
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def identify_case_specific_dna_variants(
        case_bam_file: str,
        case_bam_bai_file: str,
        control_bam_files: List[str],
        control_bam_bai_files: List[str],
        output_tsv_file: str,
        chromosomes: Optional[List[str]] = None,
        min_reads: int = CALL_DNA_VARS_MIN_READS,
        min_mapping_quality: int = CALL_DNA_VARS_MIN_MAPPING_QUALITY,
        min_average_base_quality: float = CALL_DNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        min_size_proportion: float = CALL_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance_tau: int = CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_intrachromosomal_distance: int = CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_interchromosomal_distance: int = CALL_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        apply_infinite_sites_assumption: bool = True,
        num_threads: int = CALL_DNA_VARS_NUM_THREADS,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
):
    """
    Identify case-specific DNA variants.

    Args:
        case_bam_file                       :   Case BAM file.
        case_bam_bai_file                   :   Case BAM.BAI file.
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
    if chromosomes is None:
        chromosomes = get_chromosomes(bam_file=case_bam_file)
    assert len(chromosomes) > 0
    assert(len(control_bam_files) > 0, 'control_bam_files cannot be empty.')
    assert(len(control_bam_files) == len(control_bam_bai_files), 'len(control_bam_files) must be equal to len(control_bam_bai_files).')
    df_variants = exactolibrs.identify_case_specific_dna_variants(
        case_bam_file=case_bam_file,
        case_bam_bai_file=case_bam_bai_file,
        control_bam_files=control_bam_files,
        control_bam_bai_files=control_bam_bai_files,
        output_tsv_file=output_tsv_file,
        min_reads=min_reads,
        min_mapping_quality=min_mapping_quality,
        min_average_base_quality=min_average_base_quality,
        min_size_proportion=min_size_proportion,
        max_ins_norm_edit_distance=max_ins_norm_edit_distance,
        max_intrachromosomal_distance_tau=max_intrachromosomal_distance_tau,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_interchromosomal_distance=max_interchromosomal_distance,
        apply_infinite_sites_assumption=apply_infinite_sites_assumption,
        num_threads=num_threads,
        chromosomes=chromosomes,
        temp_dir=temp_dir,
        output_type=str(output_type)
    )
    return df_variants.to_pandas()


def identify_dna_variants(
        bam_file: str,
        bam_bai_file: str,
        output_tsv_file: str,
        chromosomes: Optional[List[str]] = None,
        min_reads: int = CALL_DNA_VARS_MIN_READS,
        min_mapping_quality: int = CALL_DNA_VARS_MIN_MAPPING_QUALITY,
        min_average_base_quality: float = CALL_DNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        min_size_proportion: float = CALL_DNA_VARS_MIN_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        max_intrachromosomal_distance_tau: int = CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        max_intrachromosomal_distance: int = CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        max_interchromosomal_distance: int = CALL_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        num_threads: int = CALL_DNA_VARS_NUM_THREADS,
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
    if chromosomes is None:
        chromosomes = get_chromosomes(bam_file=bam_file)
    assert len(chromosomes) > 0
    df_variants = exactolibrs.identify_dna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        output_tsv_file=output_tsv_file,
        min_reads=min_reads,
        min_mapping_quality=min_mapping_quality,
        min_average_base_quality=min_average_base_quality,
        min_size_proportion=min_size_proportion,
        max_ins_norm_edit_distance=max_ins_norm_edit_distance,
        max_intrachromosomal_distance_tau=max_intrachromosomal_distance_tau,
        max_intrachromosomal_distance=max_intrachromosomal_distance,
        max_interchromosomal_distance=max_interchromosomal_distance,
        num_threads=num_threads,
        chromosomes=chromosomes,
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
        gene_types: List[str],
        gene_levels: List[int],
        transcript_types: List[str],
        transcript_levels: List[int],
        output_dir: str,
        output_prefix: str,
        reference_transcript_scoring_method: ReferenceTranscriptScoringMethod = ReferenceTranscriptScoringMethod(CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SCORING_METHOD),
        reference_transcript_selection_strategy: ReferenceTranscriptSelectionStrategy = ReferenceTranscriptSelectionStrategy(CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SELECTION_STRATEGY),
        reference_transcript_top_k: int = CALL_RNA_VARS_REFERENCE_TRANSCRIPT_TOP_K,
        reference_transcript_threshold: float = CALL_RNA_VARS_REFERENCE_TRANSCRIPT_THRESHOLD,
        min_mapping_quality: int = CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        min_average_base_quality: float = CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        num_threads: int = CALL_RNA_VARS_NUM_THREADS,
        temp_dir: str = "",
        output_type: OutputType = OutputType.FILE
) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame]:
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
        Pandas DataFrame of exons,
        Pandas DataFrame of read filter status,
        Pandas DataFrame of read names and transcript IDs,
        Pandas DataFrame of reference transcript matches,
        Pandas DataFrame of splice junction,
        Pandas DataFrame of transcripts,
        Pandas DataFrame of variants
    """
    df_exons,df_read_filter_status,df_read_names,df_matched_reference_transcripts,df_introns,df_transcripts,df_variant_calls = exactolibrs.identify_rna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        reference_gene_annotation_file=reference_gene_annotation_file,
        reference_gene_annotation_source=str(reference_gene_annotation_source),
        reference_gene_annotation_assembly=str(reference_gene_annotation_assembly),
        reference_gene_annotation_version=str(reference_gene_annotation_version),
        gene_types=gene_types,
        gene_levels=gene_levels,
        transcript_types=transcript_types,
        transcript_levels=transcript_levels,
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
        output_type=str(output_type)
    )
    return (df_exons.to_pandas(),
            df_read_filter_status.to_pandas(),
            df_read_names.to_pandas(),
            df_matched_reference_transcripts.to_pandas(),
            df_introns.to_pandas(),
            df_transcripts.to_pandas(),
            df_variant_calls.to_pandas())


def translate_fasta_file(
        fasta_file: str,
        temp_dir: str,
        strategy: TranslationStrategy = TranslationStrategy(TRANSLATE_STRATEGY),
        num_threads: int = TRANSLATE_NUM_THREADS
) -> pd.DataFrame:
    """
    Translate a long-read RNA-seq FASTQ file into peptide sequences.

    Args:
        fasta_file      :   FASTA file.
        strategy        :   Translation strategy.
        num_threads     :   Number of threads.

    Returns:
        Pandas DataFrame with the following columns:
        'peptide_id',
        'peptide_sequence'
        'rna_id'
        'rna_sequence'
        'orf_start'
        'orf_end'
    """
    ipc_file = exactolibrs.translate_rna_fasta_file(
        fasta_file=fasta_file,
        strategy=str(strategy),
        num_threads=num_threads,
        temp_dir=temp_dir
    )
    gc.collect()
    df = pl.read_ipc(ipc_file).to_pandas()
    if len(df) == 0:
        df = pd.DataFrame({
            'peptide_id': [],
            'peptide_sequence': [],
            'rna_id': [],
            'rna_sequence': [],
            'orf_start': [],
            'orf_end': []
        })
    return df


def translate_fastq_file(
        fastq_file: str,
        temp_dir: str,
        strategy: TranslationStrategy = TranslationStrategy(TRANSLATE_STRATEGY),
        num_threads: int = TRANSLATE_NUM_THREADS
) -> pd.DataFrame:
    """
    Translate a long-read RNA-seq FASTQ file into peptide sequences.

    Args:
        fastq_file      :   FASTA file.
        strategy        :   Translation strategy.
        num_threads     :   Number of threads.

    Returns:
        Pandas DataFrame with the following columns:
        'peptide_id',
        'peptide_sequence'
        'rna_id'
        'rna_sequence'
        'orf_start'
        'orf_end'
    """
    ipc_file = exactolibrs.translate_rna_fastq_file(
        fastq_file=fastq_file,
        strategy=str(strategy),
        num_threads=num_threads,
        temp_dir=temp_dir
    )
    gc.collect()
    df = pl.read_ipc(ipc_file).to_pandas()
    if len(df) == 0:
        df = pd.DataFrame({
            'peptide_id': [],
            'peptide_sequence': [],
            'rna_id': [],
            'rna_sequence': [],
            'orf_start': [],
            'orf_end': []
        })
    return df


def translate(rna_sequence: str, strategy: TRANSLATE_STRATEGY) -> List[Tuple[str,int,int]]:
    """
    Translate a RNA sequence.

    Args:
        rna_sequence    :   RNA sequence.
        strategy        :   Translation strategy.

    Returns:
        Peptide sequence, ORF start, ORF end
    """
    translations = exactolibrs.translate_rna_sequence(
        rna_sequence=rna_sequence,
        strategy=strategy
    )
    return translations


# def build_variation_graph(
#         variants_tsv_file: str,
#         fasta_file: str,
#         output_fasta_file: str
# ):
#     exactolibrs.build_variation_graph(variants_tsv_file, fasta_file, output_fasta_file)
#
#
# def diff_kmers(
#         query_sequences: Dict[str,str],
#         reference_sequences: Dict[str,str],
#         min_k: int,
#         max_k: int
# ) -> pd.DataFrame:
#     # Step 1. Identify all k-mers in the reference sequences
#     reference_kmers = set()
#     for sequence in reference_sequences.values():
#         for k in range(min_k, max_k + 1):
#             kmers = get_kmers(sequence=sequence, k=k)
#             reference_kmers.update(kmers)
#
#     # Step 2. Identify k-mers unique to the query sequences
#     data = {
#         'peptide_id': [],
#         'kmer': [],
#         'k': []
#     }
#     for sequence_name,sequence in query_sequences.items():
#         for k in range(min_k, max_k + 1):
#             kmers = get_kmers(sequence=sequence, k=k)
#             for kmer in kmers:
#                 if kmer not in reference_kmers:
#                     data['peptide_id'].append(sequence_name)
#                     data['kmer'].append(kmer)
#                     data['k'].append(k)
#
#     return pd.DataFrame(data)
#
#

#
#

#
# def identify_peptide_variants(
#         fasta_file: str,
#         rna_bam_file: str,
#         rna_bam_bai_file: str,
#         reference_fasta_file: str,
#         translations_tsv_file: str,
#         rna_variants_tsv_file: str,
#         dna_variants_tsv_file: str,
#         exclude_bed_file: str,
#         min_reads: int,
#         k: int,
#         num_threads: int,
#         dna_variant_padding: int,
#         output_tsv_file: str,
#         output_fasta_file: str,
#         gzip: bool
# ):
#     """
#     Identify peptide variants.
#
#     Args:
#         fasta_file             :   Input peptides FASTA file.
#         reference_fasta_file   :   Reference peptides FASTA file.
#         k                      :   K-mer k size.
#         num_threads            :   Number of threads.
#
#     Returns:
#         Pandas DataFrame with the following columns:
#         'peptide_id'
#         'peptide_sequence'
#         'kmer_mask'
#     """
#     pass
#     # exactolibrs.identify_peptide_variants(
#     #     fasta_file=fasta_file,
#     #     rna_bam_file=rna_bam_file,
#     #     rna_bam_bai_file=rna_bam_bai_file,
#     #     reference_fasta_file=reference_fasta_file,
#     #     translations_tsv_file=translations_tsv_file,
#     #     rna_variants_tsv_file=rna_variants_tsv_file,
#     #     dna_variants_tsv_file=dna_variants_tsv_file,
#     #     exclude_bed_file=exclude_bed_file,
#     #     min_reads=min_reads,
#     #     k=k,
#     #     num_threads=num_threads,
#     #     dna_variant_padding=dna_variant_padding,
#     #     output_tsv_file=output_tsv_file,
#     #     output_fasta_file=output_fasta_file,
#     #     gzip=gzip,
#     #     output_type='file'
#     # )
#
#
# def integrate_dna_rna_variants(
#         annotated_dna_variant_callset_tsv_file: str,
#         rna_variant_callset_tsv_file: str,
#         reference_gene_annotation_file: str,
#         reference_gene_annotation_source: GeneAnnotationSource,
#         output_tsv_file: str,
#         max_exon_offset: int = INTEGRATE_VARS_MAX_EXON_OFFSET,
#         max_transcript_boundary_offset: int = INTEGRATE_VARS_MAX_TRANSCRIPT_BOUNDARY_OFFSET,
#         max_intergenic_distance: int = INTEGRATE_VARS_MAX_INTERGENIC_DISTANCE,
#         num_threads: int = CALL_RNA_VARS_NUM_THREADS,
#         temp_dir: str = "",
#         output_type: OutputType = OutputType.FILE
# ) -> pd.DataFrame:
#     """
#     Integrate DNA and RNA variants.
#
#     Args:
#         annotated_dna_variant_callset_tsv_file      :   Annotated DNA variant callset TSV file.
#         rna_variant_callset_tsv_file                :   RNA variant callset TSV file.
#         reference_gene_annotation_file              :   Reference gene annotation TSV file.
#         reference_gene_annotation_source            :   Reference gene annotation TSV file source.
#         output_tsv_file                             :   Output TSV file.
#         max_exon_offset                             :   Maximum exon offset.
#         max_transcript_boundary_offset              :   Maximum transcript boundary offset.
#         max_intergenic_distance                     :   Maximum intergenic distance.
#         num_threads                                 :   Number of threads.
#         temp_dir                                    :   Temp directory (default: TMPDIR).
#         output_type                                 :   Output type ('file' or 'dataframe').
#
#     Returns:
#         Pandas DataFrame with the following columns:
#             'rna_variant_call_id,
#             'dna_variant_call_id',
#             'distance',
#             'rna_variant_position',
#             'dna_variant_position'
#     """
#     df_integration = exactolibrs.integrate_dna_rna_variants(
#         annotated_dna_variant_callset_tsv_file=annotated_dna_variant_callset_tsv_file,
#         rna_variant_callset_tsv_file=rna_variant_callset_tsv_file,
#         reference_gene_annotation_file=reference_gene_annotation_file,
#         reference_gene_annotation_source=str(reference_gene_annotation_source),
#         output_tsv_file=output_tsv_file,
#         max_exon_offset=max_exon_offset,
#         max_transcript_boundary_offset=max_transcript_boundary_offset,
#         max_intergenic_distance=max_intergenic_distance,
#         num_threads=num_threads,
#         temp_dir=temp_dir,
#         output_type=str(output_type)
#     )
#     return df_integration.to_pandas()
#
#
# def remove_unspliced_rnas(
#         bam_file: str,
#         fasta_file: str,
#         reference_gene_annotation_file: str,
#         reference_gene_annotation_source: GeneAnnotationSource,
#         output_bam_file: str,
#         output_fasta_file: str,
#         output_tsv_file: str,
#         num_threads: int,
#         min_mapping_quality: int
# ):
#     if reference_gene_annotation_source == GeneAnnotationSource.GENCODE:
#         # Get single-exon transcripts
#         gencode = Gencode(gtf_file=gencode_gtf_file, version='', species='')
#         exon_counts = gencode.df_exons.groupby('transcript_id')['exon_id'].nunique()
#         single_exon_transcript_ids = exon_counts[exon_counts == 1].index
#         df_transcripts = gencode.df_transcripts[
#             gencode.df_transcripts['transcript_id'].isin(single_exon_transcript_ids)
#         ]
#         return df_transcripts
#     else:
#         raise Exception("Unsupported reference gene annotation source: %s" % reference_gene_annotation_source)
#
#
#
