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


import pandas as pd
import pysam
from dataclasses import field
from typing import List
from .alignment_map import AlignmentMap
from .annotation_db import AnnotationDb
from .constants import VariantCallingMethods, NucleicAcidTypes
from .default_parameters import *
from .fasta import Fasta
from .gene_set import GeneSet
from .genomic_ranges_list import  GenomicRangesList
from .logging import get_logger
from .variant_filter import VariantFilter
from .variants_list import VariantsList
from .vcf import Vcf


logger = get_logger(__name__)


def run_exacto_vcf_to_tsv(
        vcf_file: str,
        source_id: str,
        variant_calling_method: str,
        sequencing_platform: str
) -> VariantsList:
    """
    Convert a VCF file to a VariantsList.

    Parameters
    ----------
    vcf_file                :   VCF file.
    source_id               :   Source ID (e.g. patient ID or cell line sample ID).
    variant_calling_method  :   Variant calling method.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    variants_list           :   VariantsList object.
    """
    df_vcf = Vcf.read_vcf_file(vcf_file=vcf_file)
    if variant_calling_method == VariantCallingMethods.CUTESV:
        variants_list = Vcf.parse_cutesv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        variants_list = Vcf.parse_deepvariant_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        variants_list = Vcf.parse_gatk4_mutect2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.PBSV:
        variants_list = Vcf.parse_pbsv_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.SNIFFLES2:
        variants_list = Vcf.parse_sniffles2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.STRELKA2:
        variants_list = Vcf.parse_strelka2_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    elif variant_calling_method == VariantCallingMethods.SVIM:
        variants_list = Vcf.parse_svim_callset(
            df_vcf=df_vcf,
            sequencing_platform=sequencing_platform,
            source_id=source_id
        )
    else:
        raise Exception('Unsupported variant calling method: %s' % variant_calling_method)
    return variants_list


def run_exacto_merge_variant_calls(
        variants_lists: List[VariantsList],
        max_neighbor_distance: int = MERGE_MAX_NEIGHBOR_DISTANCE,
) -> VariantsList:
    """
    Merges VariantsList objects into one.

    Parameters
    ----------
    variants_lists                  :   List of VariantsList objects.
    max_neighbor_distance           :   Maximum neighbor distance.

    Returns
    -------
    variants_list                   :   VariantsList object.
    """
    variants_list = VariantsList.merge(
        variants_lists=variants_lists,
        max_neighbor_distance=max_neighbor_distance
    )
    return variants_list


def run_exacto_filter_variants(
        variants_list: VariantsList,
        variant_filters: List[VariantFilter] = field(default_factory=list),
        excluded_variants_list: VariantsList = None,
        excluded_regions_list: GenomicRangesList = None,
        excluded_variants_padding: int = FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING,
        excluded_regions_padding: int = FILTER_VARIANTS_EXCLUDED_REGION_PADDING,
        num_processes: int = FILTER_VARIANTS_NUM_PROCESSES
) -> VariantsList:
    """
    Filters a VariantsList object.

    Parameters
    ----------
    variants_list               :   VariantsList object.
    variant_filters             :   List of VariantFilter objects.
    excluded_variants_list      :   VariantsList object of variants to exclude.
    excluded_regions_list       :   GenomicRangesList object of regions to exclude.
    excluded_variants_padding   :   Number of bases to pad each variant's positions 1 and 2.
    excluded_regions_padding    :   Number of bases to pad each region to exclude.
    num_processes               :   Number of processes.

    Returns
    -------
    variants_list               :   VariantsList object.
    """
    logger.info('%i variants in the original list before any filtering.' % variants_list.size)
    logger.info('%i variant calls in the original list before any filtering.' % len(variants_list.variant_call_ids))

    # Step 1. Filter out variants based on VariantFilter
    if len(variant_filters) > 0:
        variants = variants_list.filter(
            variant_filters=variant_filters,
            num_processes=num_processes,
        )
        variant_ids = [variant.id for variant in variants]
        for variant in variants_list.variants:
            if variant.id not in variant_ids:
                variants_list.remove_by_id(id=variant.id)
        logger.info('%i variants remain after applying variant filters.' % variants_list.size)
        logger.info('%i variant calls remain after applying variant filters.' % len(variants_list.variant_call_ids))

    # Step 2. Filter out variants near the excluded variants
    if excluded_variants_list is not None:
        nearby_variants = variants_list.find_nearby_variants(
            variants=excluded_variants_list.variants,
            padding=excluded_variants_padding,
            num_processes=num_processes
        )
        for nearby_variant in nearby_variants:
            variants_list.remove_by_id(id=nearby_variant.id)
        logger.info('%i variants remain after removing variants near excluded variants.' % variants_list.size)
        logger.info('%i variant calls remain after removing variants near excluded variants.' % len(variants_list.variant_call_ids))

    # Step 3. Filter out variants near the excluded regions
    if excluded_regions_list is not None:
        nearby_variants = variants_list.filter_regions(
            genomic_ranges_list=excluded_regions_list,
            padding=excluded_regions_padding,
            num_processes=num_processes
        )
        for nearby_variant in nearby_variants:
            variants_list.remove_by_id(id=nearby_variant.id)
        logger.info('%i variants remain after removing variants near excluded regions.' % variants_list.size)
        logger.info('%i variant calls remain after removing variants near excluded regions.' % len(variants_list.variant_call_ids))

    return variants_list


def run_exacto_annotate_variants_list(
        variants_list: VariantsList,
        annotation_db: AnnotationDb
) -> VariantsList:
    """
    Annotates a variants list and returns the annotated variants list.

    Parameters
    ----------
    variants_list       :   VariantsList object.
    annotation_db       :   AnnotationDb object.

    Returns
    -------
    variants_list       :   VariantsList object.
    """
    return annotation_db.annotate_variants_list(variants_list=variants_list)


def run_exacto_simulate_variants(
        genome_fasta: pysam.FastaFile,
        gene_set: GeneSet,
        df_target_regions: pd.DataFrame,
        num_snv: int = SIMULATE_VARIANTS_NUM_SNV,
        num_insertion: int = SIMULATE_VARIANTS_NUM_INSERTION,
        num_deletion: int = SIMULATE_VARIANTS_NUM_DELETION,
        insertion_size_mean: int = SIMULATE_VARIANTS_INSERTION_SIZE_MEAN,
        insertion_size_stdev: int = SIMULATE_VARIANTS_INSERTION_SIZE_STDEV,
        deletion_size_mean: int = SIMULATE_VARIANTS_DELETION_MEAN,
        deletion_size_stdev: int = SIMULATE_VARIANTS_DELETION_STDEV,
        enforce_infinite_sites_model: bool = SIMULATE_VARIANTS_ENFORCE_INFINITE_SITES_MODEL
    ) -> pd.DataFrame:
    """
    Simulates DNA and RNA variants.

    Parameters
    ----------
    genome_fasta                            :   pysam.FastaFile object of reference genome.
    gene_set                                :   An instance of 'GeneSet' class.
    num_snv                                 :   Number of SNVs to simulate.
    num_insertion                           :   Number of insertions to simulate.
    num_deletion                            :   Number of deletions to simulate.
    insertion_size_mean                     :   Mean value of insertion size.
    insertion_size_stdev                    :   Standard deviation of insertion size.
    deletion_size_mean                      :   Mean value of deletion size.
    deletion_size_stdev                     :   Standard deviation of deletion size.
    enforce_infinite_sites_model            :   If true, the simulation enforces the infinite sites model.

    Returns
    -------
    df_variants                             :   DataFrame of variants.
    """
    a = 1
    # df_rna_variants, variant_transcript_sequences = simulate_rna_variants(**locals())
    # return df_rna_variants, variant_transcript_sequences


def run_exacto_call_rna_variants(
        bam: pysam.AlignmentFile,
        num_processes: int
    ) -> pd.DataFrame:
    """
    Call RNA variants in a BAM file.

    Parameters
    ----------
    bam                 :   Pysam AlignmentFile object of a BAM file.
    num_processes       :   Number of processes.

    Returns
    -------
    variants_list       :   An instance of 'VariantsList' class.
    """
    alignment_map = AlignmentMap(bam=bam, nucleic_acid_type=NucleicAcidTypes.RNA)
    return alignment_map.call_variants(num_processes=num_processes)
    # bam_filename = str(bam_file.filename.decode())
    # variant_callset = exactors.identify_rna_variants(bam_filename, num_cores)
    # df_variants = pd.DataFrame({
    #     'chromosome': variant_callset.chromosomes,
    #     'position': variant_callset.positions,
    #     'read_id': variant_callset.read_ids,
    #     'variant_type': variant_callset.variant_types,
    #     'reference_allele': variant_callset.reference_alleles,
    #     'alternate_allele': variant_callset.alternate_alleles,
    #     'sequence': variant_callset.sequences,
    #     'variant_size': variant_callset.variant_sizes
    # })
    # return df_variants
#
#
# def run_exacto_simulate_rna_variants(
#         genome_fasta: pysam.FastaFile,
#         df_genes: pd.DataFrame,
#         df_transcripts: pd.DataFrame,
#         df_exons: pd.DataFrame,
#         df_target_regions: pd.DataFrame,
#         df_herv_regions: pd.DataFrame,
#         num_snv: int = SIMULATE_RNA_VARIANTS_NUM_SNV,
#         num_insertion: int = SIMULATE_RNA_VARIANTS_NUM_INSERTION,
#         num_deletion: int = SIMULATE_RNA_VARIANTS_NUM_DELETION,
#         num_fusion: int = SIMULATE_RNA_VARIANTS_NUM_FUSION,
#         num_inversion: int = SIMULATE_RNA_VARIANTS_NUM_INVERSION,
#         num_herv: int = SIMULATE_RNA_VARIANTS_NUM_HERV,
#         insertion_size_mean: int = SIMULATE_RNA_VARIANTS_INSERTION_SIZE_MEAN,
#         insertion_size_stdev: int = SIMULATE_RNA_VARIANTS_INSERTION_SIZE_STDEV,
#         deletion_size_mean: int = SIMULATE_RNA_VARIANTS_DELETION_MEAN,
#         deletion_size_stdev: int = SIMULATE_RNA_VARIANTS_DELETION_STDEV,
#         herv_solo_ltr_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_SOLO_LTR,
#         herv_truncated_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_TRUNCATED,
#         herv_chimeric_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_CHIMERIC,
#         herv_chimeric_max_neighboring_distance: int = SIMULATE_RNA_VARIANTS_HERV_CHIMERIC_MAX_NEIGHBORING_DISTANCE,
#         herv_full_length_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_FULL_LENGTH,
#         infinite_sites_assumption: bool = SIMULATE_RNA_VARIANTS_INFINITE_SITES_ASSUMPTION
#     ) -> Tuple[pd.DataFrame, List]:
#     """
#     Simulates RNA variants.
#
#     Parameters
#     ----------
#     genome_fasta                            :   pysam.FastaFile object of reference genome.
#     df_transcripts                          :   DataFrame of transcripts.
#     df_exons                                :   DataFrame of exons.
#     df_target_regions                       :   DataFrame of regions to simulate RNA variants.
#     df_herv_regions                         :   HERV regions.
#     num_snv                                 :   Number of SNVs to simulate.
#     num_insertion                           :   Number of insertions to simulate.
#     num_deletion                            :   Number of deletions to simulate.
#     num_fusion                              :   Number of fusions to simulate.
#     num_inversion                           :   Number of inversions to simulate.
#     num_herv                                :   Number of HERVs to simulate.
#     insertion_size_mean                     :   Mean value of insertion size.
#     insertion_size_stdev                    :   Standard deviation of insertion size.
#     deletion_size_mean                      :   Mean value of deletion size.
#     deletion_size_stdev                     :   Standard deviation of deletion size.
#     herv_solo_ltr_proportion                :   Proportion of expressed HERVs that only have solo LTR sequences.
#     herv_truncated_proportion               :   Proportion of HERVs that are truncated.
#     herv_chimeric_proportion                :   Proportion of HERVs that are chimeric (concatenation of neighboring HERVs).
#     herv_chimeric_max_neighboring_distance  :   Maximum distance for two HERVs to be considered for simulation of a
#                                                 chimeric HERV.
#     herv_full_length_proportion             :   Proportion of HERVs that are full-lengths.
#     infinite_sites_assumption               :   If true, the simulation enforces the infinite sites assumption.
#
#     Returns
#     -------
#     df_rna_variants                         :   DataFrame of RNA variants
#     variant_transcript_sequences            :   List of variant transcript sequences
#     """
#     df_rna_variants, variant_transcript_sequences = simulate_rna_variants(**locals())
#     return df_rna_variants, variant_transcript_sequences
#
#
def run_exacto_simulate_reads(
        fasta: Fasta,
        output_fastq_gz_file: str,
        num_gigabases: float,
        read_length_mean: float,
        read_length_stdev: float,
        base_quality_mean: float,
        base_quality_stdev: float
    ) -> None:
    """
    Simulates sequencing reads.

    Parameters
    ----------
    sequences               :   List of instances of the class Sequence.
    output_fastq_gz_file    :   Output .fastq.gz file.
    num_gigabases           :   Number of gigabases to sequence.
    read_length_mean        :   Mean value of read length.
    read_length_stdev       :   Standard deviation of read length.
    base_quality_mean       :   Mean value of base quality.
    base_quality_stdev      :   Standard deviation of base quality.

    Returns
    -------
    reads               :   List of instances of the class Read.
    """
    pass
    # return simulate_single_end_reads(sequences=sequences,
    #                                  output_fastq_gz_file=output_fastq_gz_file,
    #                                  num_bases=num_gigabases * 10e9,
    #                                  read_length_mean=read_length_mean,
    #                                  read_length_stdev=read_length_stdev,
    #                                  base_quality_mean=base_quality_mean,
    #                                  base_quality_stdev=base_quality_stdev)
