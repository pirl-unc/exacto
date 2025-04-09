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
from typing import List, Optional, Tuple
from .default import *
from .logging import get_logger
from .utilities import get_chromosomes


logger = get_logger(__name__)


def identify_dna_variants(
        bam_file: str,
        bam_bai_file: str,
        output_tsv_file: str,
        chromosomes: Optional[List[str]] = None,
        gzip: bool = True,
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
        output_type: str = "file"
) -> pd.DataFrame:
    """
    Identify DNA variants.

    Args:
        bam_file                            :   BAM file.
        bam_bai_file                        :   BAM.BAI file.
        output_tsv_file                     :   Output TSV file.
        chromosomes                         :   A list of chromosomes to scan for the presence of DNA variants (default: None).
                                                If left unspecified (i.e. None), all chromosomes in the BAM file will be considered.
        gzip                                :   If True, gzip the output TSV file.
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
    assert output_type in ["file", "dataframe"]
    df_variants = exactolibrs.identify_dna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        output_tsv_file=output_tsv_file,
        gzip=gzip,
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
        output_type=output_type
    )
    return df_variants.to_pandas()


def identify_case_specific_dna_variants(
        case_bam_file: str,
        case_bam_bai_file: str,
        control_bam_files: List[str],
        control_bam_bai_files: List[str],
        output_tsv_file: str,
        chromosomes: Optional[List[str]] = None,
        gzip: bool = True,
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
        output_type: str = "file"
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
        gzip                                :   If True, gzip the output TSV file.
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
    assert output_type in ["file", "dataframe"]
    df_variants = exactolibrs.identify_case_specific_dna_variants(
        case_bam_file=case_bam_file,
        case_bam_bai_file=case_bam_bai_file,
        control_bam_files=control_bam_files,
        control_bam_bai_files=control_bam_bai_files,
        output_tsv_file=output_tsv_file,
        gzip=gzip,
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
        output_type=output_type
    )
    return df_variants.to_pandas()


def identify_rna_variants(
        bam_file: str,
        bam_bai_file: str,
        reference_genome_fasta_file: str,
        gene_annotation_file: str,
        gene_annotation_source: str,
        output_exons_tsv_file: str,
        output_sj_tsv_file: str,
        output_variants_tsv_file: str,
        gzip: bool = True,
        min_mapping_quality: int = CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        min_average_base_quality: float = CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        num_threads: int = CALL_RNA_VARS_NUM_THREADS,
        temp_dir: str = "",
        output_type: str = "file"
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
        If output_type is 'dataframe',
        then exons Pandas DataFrame, splice junction Pandas DataFrame, variants Pandas DataFrame.
    """
    df_exons,df_splice_junctions,df_variants = exactolibrs.identify_rna_variants(
        bam_file=bam_file,
        bam_bai_file=bam_bai_file,
        reference_genome_fasta_file=reference_genome_fasta_file,
        gene_annotation_file=gene_annotation_file,
        gene_annotation_source=gene_annotation_source,
        output_exons_tsv_file=output_exons_tsv_file,
        output_splice_junctions_tsv_file=output_sj_tsv_file,
        output_variants_tsv_file=output_variants_tsv_file,
        gzip=gzip,
        min_mapping_quality=min_mapping_quality,
        min_average_base_quality=min_average_base_quality,
        num_threads=num_threads,
        temp_dir=temp_dir,
        output_type=output_type
    )
    return df_exons.to_pandas(), df_splice_junctions.to_pandas(), df_variants.to_pandas()


def identify_peptide_variants(
        fasta_file: str,
        rna_bam_file: str,
        rna_bam_bai_file: str,
        reference_fasta_file: str,
        translations_tsv_file: str,
        rna_variants_tsv_file: str,
        dna_variants_tsv_file: str,
        exclude_bed_file: str,
        min_reads: int,
        k: int,
        num_threads: int,
        dna_variant_padding: int,
        output_tsv_file: str,
        output_fasta_file: str,
        gzip: bool
):
    """
    Identify peptide variants.

    Args:
        fasta_file             :   Input peptides FASTA file.
        reference_fasta_file   :   Reference peptides FASTA file.
        k                      :   K-mer k size.
        num_threads            :   Number of threads.

    Returns:
        Pandas DataFrame with the following columns:
        'peptide_id'
        'peptide_sequence'
        'kmer_mask'
    """
    pass
    # exactolibrs.identify_peptide_variants(
    #     fasta_file=fasta_file,
    #     rna_bam_file=rna_bam_file,
    #     rna_bam_bai_file=rna_bam_bai_file,
    #     reference_fasta_file=reference_fasta_file,
    #     translations_tsv_file=translations_tsv_file,
    #     rna_variants_tsv_file=rna_variants_tsv_file,
    #     dna_variants_tsv_file=dna_variants_tsv_file,
    #     exclude_bed_file=exclude_bed_file,
    #     min_reads=min_reads,
    #     k=k,
    #     num_threads=num_threads,
    #     dna_variant_padding=dna_variant_padding,
    #     output_tsv_file=output_tsv_file,
    #     output_fasta_file=output_fasta_file,
    #     gzip=gzip,
    #     output_type='file'
    # )


def translate(rna_sequence: str, strategy: TRANSLATE_STRATEGY):
    """
    Translate a RNA sequence.

    Args:
        rna_sequence    :   RNA sequence.
        strategy        :   Translation strategy.

    Returns:
        Peptide sequence, ORF start, ORF end
    """
    sequence, orf_start, orf_end = exactolibrs.translate_rna_sequence(rna_sequence=rna_sequence, strategy=strategy)
    return sequence, orf_start, orf_end


def translate_fasta_file(
        fasta_file: str,
        temp_dir: str,
        strategy: TRANSLATE_STRATEGY,
        num_threads: int = TRANSLATE_NUM_THREADS,
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
        strategy=strategy,
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
        strategy: TRANSLATE_STRATEGY,
        num_threads: int = TRANSLATE_NUM_THREADS,
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
        strategy=strategy,
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
