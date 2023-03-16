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
    variants_lists = VariantsList.merge(
        variants_lists=variants_lists,
        enforce_variant_type_matching=enforce_variant_type_matching,
        max_neighbor_distance=max_neighbor_distance
    )
    return variants_lists

#
# def run_exacto_refine(
#         variants_list: VariantsList,
#         df_variants_to_exclude: pd.DataFrame,
#         df_gapped_regions: pd.DataFrame,
#         queries: List[str] = [],
#         gapped_regions_padding: int = GENOME_GAPPED_REGIONS_PADDING,
#         exclude_variants_padding: int = EXCLUDE_SV_PADDING,
#         num_processes: int = NUM_PROCESSES_REFINE
#     ) -> pd.DataFrame:
#     """
#     Refines a set of structural variants and returns the refined set.
#
#     Parameters
#     ----------
#     df_variants                         :   DataFrame of structural variants.
#                                             Expected columns:
#                                             'chr_1'
#                                             'pos_1'
#                                             'chr_2'
#                                             'pos_2'
#                                             'variant_type'
#     df_variants_to_exclude              :   DataFrame of variants to exclude.
#     df_gapped_regions                   :   DataFrame of gapped regions in the genome.
#     gapped_regions_padding              :   Number of bases to pad gapped regions'
#                                             start and end positions.
#     exclude_variants_padding            :   Number of bases to pad breakpoints of
#                                             structural variants in df_structural_variants_to_exclude
#     num_processes                       :   Number of processes.
#
#     Returns
#     -------
#     DataFrame of refined variants.
#     """
#     logger.info('%i variants before refinement.' % len(df_variants))
#
#     # Filter out variants based on queries
#     df_variants = refine_variants(df_variants=df_variants, queries=queries)
#
#     # Filter out variants near the gapped regions.
#     if df_gapped_regions is not None:
#         df_variants = remove_variants_near_gapped_regions(
#             df_variants=df_variants,
#             df_gapped_regions=df_gapped_regions,
#             gapped_regions_padding=gapped_regions_padding,
#             num_processes=num_processes
#         )
#         logger.info(
#             '%i variants after filtering out variants near the gapped regions.'
#             % len(df_variants)
#         )
#
#     # Filter out excluded variants.
#     if df_variants_to_exclude is not None:
#         df_variants = remove_variants(
#             df_variants=df_variants,
#             df_variants_to_exclude=df_variants_to_exclude,
#             exclude_variants_padding=exclude_variants_padding,
#             num_processes=num_processes
#         )
#         logger.info('%i variants after filtering out excluded variants.' % len(df_variants))
#
#     logger.info('%i variants after refinement.' % len(df_variants))
#     return df_variants
#
#
# def run_exacto_annotate(
#         df_variants: pd.DataFrame,
#         annotation_source: str,
#         df_gencode_genes: pd.DataFrame,
#         df_gencode_exons: pd.DataFrame,
#         ensembl_release: int = -1,
#         ensembl_species: str = ''
#     ) -> pd.DataFrame:
#     """
#     Annotates a set of variants and returns the annotated set.
#
#     Parameters
#     ----------
#     df_variants         :   DataFrame of variants.
#                             Expected columns:
#                             'chr_1'
#                             'pos_1'
#                             'chr_2'
#                             'pos_2'
#                             'sv_type' (DEL, INS, INV, DUP, BND or TRA)
#     annotation_source   :   Annotation source ('ensembl' or 'gencode').
#     df_gencode_genes    :   DataFrame of GENCODE genes.
#                             Specify this if 'annotation_source' is 'gencode'.
#                             Expected columns:
#                             'gene_id'
#                             'gene_name'
#                             'gene_type'
#                             'gene_chrom'
#                             'gene_start'
#                             'gene_end'
#                             'gene_strand'
#                             'level'
#                             'transcripts_count'
#     df_gencode_exons    :   DataFrame of GENCODE exons.
#                             Specify this if 'annotation_source' is 'gencode'.
#                             Expected columns:
#                             'gene_id'
#                             'transcript_id'
#                             'exon_id'
#                             'exon_number'
#                             'exon_chrom'
#                             'exon_start'
#                             'exon_end'
#     ensembl_release     :   Ensembl release version.
#                             Specify this if 'annotation_source' is 'ensembl'.
#     ensembl_species     :   Ensembl species.
#
#     Returns
#     -------
#     DataFrame of annotated genomic structural variants.
#     """
#     if len(df_variants) == 0:
#         logger.warning('DataFrame is empty. Returning without annotating.')
#         return df_variants
#     if annotation_source == AnnotationSources.ENSEMBL:
#         df_structural_variants = annotate_variants_using_pyensembl(
#             df_variants=df_variants,
#             ensembl_release=ensembl_release,
#             species=ensembl_species
#         )
#     elif annotation_source == AnnotationSources.GENCODE:
#         df_structural_variants = annotate_variants_using_gencode(
#             df_variants=df_variants,
#             df_gencode_genes=df_gencode_genes,
#             df_gencode_exons=df_gencode_exons
#         )
#     else:
#         raise Exception(
#             "Invalid value for 'annotation_source': %s. Allowed 'annotation_source' values are %s "
#             % (annotation_source, ', '.join(AnnotationSources.ALL)))
#     return df_structural_variants
#
#

#
#
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
