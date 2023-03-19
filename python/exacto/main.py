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
from typing import List
from exacto import exactors
from .vcf import read_vcf_file
from .annotation import Annotation
from .ensembl import Ensembl
from .gencode import Gencode
from .variant_filter import VariantFilter
from .variants_list import VariantsList
from .logging import get_logger
from .default_parameters import *


logger = get_logger(__name__)


def run_exacto_convert(
        vcf_file: str,
        source_id: str,
        variant_calling_method: str,
        sequencing_platform: str,
        tumor_sample_id: str,
        normal_sample_id: str
    ) -> VariantsList:
    """
    Converts a VCF file to a VariantsList.

    Parameters
    ----------
    vcf_file                :   VCF file.
    source_id               :   Source ID (e.g. patient ID or cell line sample ID).
    variant_calling_method  :   Variant calling method.
    sequencing_platform     :   Sequencing platform.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID.

    Returns
    -------
    variants_list           :   An instance of the class VariantsList.
    """
    variants_list = read_vcf_file(
        vcf_file=vcf_file,
        variant_calling_method=variant_calling_method,
        sequencing_platform=sequencing_platform,
        source_id=source_id,
        tumor_sample_id=tumor_sample_id,
        normal_sample_id=normal_sample_id
    )
    return variants_list


def run_exacto_merge(
        variants_lists: List[VariantsList],
        enforce_variant_type_matching: bool = True,
        max_neighbor_distance: int = MAX_NEIGHBOR_DISTANCE,
    ) -> VariantsList:
    """
    Merges a list of variant.

    Parameters
    ----------
    variants_lists                  :   List of instances of the VariantsList class.
    enforce_variant_type_matching   :   If true, variant_type must match for two VariantCalls
                                        to be in the same Variant.
    max_neighbor_distance           :   Maximum neighbor distance.

    Returns
    -------
    variants_list                   :   An instance of the VariantsList class.
    """
    variants_list = VariantsList.merge(
        variants_lists=variants_lists,
        enforce_variant_type_matching=enforce_variant_type_matching,
        max_neighbor_distance=max_neighbor_distance
    )
    return variants_list


def run_exacto_filter(
        variants_list: VariantsList,
        df_excluded_variants: pd.DataFrame,
        df_excluded_regions: pd.DataFrame,
        variant_filters: List[VariantFilter],
        excluded_region_padding: int = EXCLUDED_REGION_PADDING,
        excluded_variant_padding: int = EXCLUDED_VARIANT_PADDING,
        enforce_variant_type_matching: bool = ENFORCE_VARIANT_TYPE_MATCHING,
        num_processes: int = NUM_PROCESSES_FILTER
    ) -> VariantsList:
    """
    Filters a variants list.

    Parameters
    ----------
    variants_list                   :   An instance of the VariantsList class.
    df_excluded_variants            :   DataFrame. Expected columns:
                                        'chr_1', 'pos_1', 'chr_2', 'pos_2'
    df_excluded_regions             :   DataFrame. Expected columns:
                                        'chrom', 'chromStart', 'chromEnd'
    variant_filters                 :   List of instances of the VariantFilter class.
    excluded_region_padding         :   Number of bases to pad each region to exclude.
    excluded_variant_padding        :   Number of bases to pad each variant's positions 1 and 2.
    enforce_variant_type_matching   :   Enforce variant type matching.
    num_processes                   :   Number of processes.

    Returns
    -------
    variants_list                   :   An instance of the VariantsList class.
    """
    logger.info('%i variants in the original list before filtering.' % len(variants_list.variant_ids))
    logger.info('%i variant calls in the original list before filtering.' % len(variants_list.variant_call_ids))

    # Step .1 Filter out variants based on variant filters
    if len(variant_filters) > 0:
        variants_list.filter(
            variant_filters=variant_filters,
            num_processes=num_processes
        )
        logger.info('%i variants remain after applying variant filters.' % len(variants_list.variant_ids))
        logger.info('%i variant calls remain after applying variant filters.' % len(variants_list.variant_call_ids))

    # Step 2. Filter out variants near the excluded regions.
    if len(df_excluded_regions) > 0:
        variants_list.filter_regions(
            df_excluded_regions=df_excluded_regions,
            excluded_regions_padding=excluded_region_padding
        )
        logger.info('%i variants remain after removing variant calls near excluded regions.' % len(variants_list.variant_ids))
        logger.info('%i variant calls remain after removing variant calls near excluded regions.' % len(variants_list.variant_call_ids))

    # Step 3. Filter out variants near the excluded variants
    if len(df_excluded_variants) > 0:
        variants_list.filter_variants(
            df_excluded_variants=df_excluded_variants,
            excluded_variant_padding=excluded_variant_padding,
            enforce_variant_type_matching=enforce_variant_type_matching
        )
        logger.info('%i variants remain after removing variant calls near excluded variants.' % len(variants_list.variant_ids))
        logger.info('%i variant calls remain after removing variant calls near excluded variants.' % len(variants_list.variant_call_ids))

    return variants_list


def run_exacto_annotate(
        variants_list: VariantsList,
        annotation: Annotation,
    ) -> pd.DataFrame:
    """
    Annotates a variants list and returns the annotated variants list.

    Parameters
    ----------
    variants_list       :   An instance of the VariantsList class.
    annotation          :   An instance of the Annotation class.

    Returns
    -------
    variants_list       :   An instance of the VariantsList class.
    """
    variants_list = annotation.annotate_variants(variants_list=variants_list)
    return variants_list


# def run_exacto_identify_rna_variants(
#         bam_file: pysam.AlignmentFile,
#         num_cores: int
#     ) -> pd.DataFrame:
#     """
#     Identifies RNA variants in a BAM file.
#
#     Parameters
#     ----------
#     bam_file    :   Pysam AlignmentFile object of a BAM file.
#     num_cores   :   Number of cores to use.
#
#     Returns
#     -------
#     """
#     bam_filename = str(bam_file.filename.decode())
#     variant_callset = exactors.identify_rna_variants(bam_filename, num_cores)
#     df_variants = pd.DataFrame({
#         'chromosome': variant_callset.chromosomes,
#         'position': variant_callset.positions,
#         'read_id': variant_callset.read_ids,
#         'variant_type': variant_callset.variant_types,
#         'reference_allele': variant_callset.reference_alleles,
#         'alternate_allele': variant_callset.alternate_alleles,
#         'sequence': variant_callset.sequences,
#         'variant_size': variant_callset.variant_sizes
#     })
#     return df_variants
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
# def run_exacto_simulate_reads(
#         sequences: List[Sequence],
#         num_gigabases: float,
#         read_length_mean: float,
#         read_length_stdev: float,
#         base_quality_mean: float,
#         base_quality_stdev: float
#     ) -> List[Read]:
#     """
#     Simulates sequencing reads.
#
#     Parameters
#     ----------
#     sequences           :   List of instances of the class Sequence.
#     num_gigabases       :   Number of gigabases to sequence.
#     read_length_mean    :   Mean value of read length.
#     read_length_stdev   :   Standard deviation of read length.
#     base_quality_mean   :   Mean value of base quality.
#     base_quality_stdev  :   Standard deviation of base quality.
#
#     Returns
#     -------
#     reads               :   List of instances of the class Read.
#     """
#     return simulate_single_end_reads(sequences=sequences,
#                                      num_bases=num_gigabases * 10e9,
#                                      read_length_mean=read_length_mean,
#                                      read_length_stdev=read_length_stdev,
#                                      base_quality_mean=base_quality_mean,
#                                      base_quality_stdev=base_quality_stdev)
