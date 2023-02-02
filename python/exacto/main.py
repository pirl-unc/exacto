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


from __future__ import print_function, division, absolute_import


import pysam
import random
import pandas as pd
from typing import Tuple, List
from .constants import *
from .default_parameters import *
from .logging import get_logger
from .variant_refinement.structural_variants import *
from .variant_refinement.small_variants import *
from .variant_annotation.ensembl import *
from .variant_annotation.gencode import *
from .variant_annotation.annovar import *
from .utilities.merging_utils import *
from .utilities.gencode_utils import *
from .simulation.rna_variants import *
from .simulation.reads import *
from exacto import exactors


logger = get_logger(__name__)


def run_exacto_refine_genomic_structural_variants(
        df_structural_variants: pd.DataFrame,
        df_structural_variants_to_exclude: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        variant_calling_method: str,
        keep_only_precise_sv: bool,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding: int = GENOME_GAPPED_REGIONS_PADDING,
        exclude_variants_padding: int = EXCLUDE_SV_PADDING) -> pd.DataFrame:
    """
    Refines a set of structural variants and returns the refined set.

    Parameters
    ----------
    df_structural_variants              :   DataFrame of structural variants.
                                            Expected columns:

    df_structural_variants_to_exclude   :   DataFrame of structural variants to exclude.
    df_gapped_regions                   :   DataFrame of gapped regions in the genome.
    variant_calling_method              :   Variant calling method.
    keep_only_precise_sv                :   If true, only structural variants with
                                            'precise' breakpoints are kept.
    keep_only_chromosomes               :   List of chromosomes to keep.
                                            Chromosomes not specified in this list
                                            will be removed.
    keep_only_filter_values             :   List of 'filter' values to keep.
                                            'filter' values not specified in this list
                                            will be removed.
    min_total_depth                     :   Minimum total depth.
    min_variant_reads_count             :   Minimum variant reads count.
    gapped_regions_padding              :   Number of bases to pad gapped regions'
                                            start and end positions.
    exclude_variants_padding            :   Number of bases to pad breakpoints of
                                            structural variants in df_structural_variants_to_exclude

    Returns
    -------
    DataFrame of refined structural variants.
    """
    logger.info('%i variants before refinement.' % len(df_structural_variants))

    if variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2:
        df_structural_variants = refine_sniffles2_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            keep_only_precise=keep_only_precise_sv,
            min_total_depth=min_total_depth,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.CUTESV:
        df_structural_variants = refine_cutesv_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            keep_only_precise=keep_only_precise_sv,
            min_total_depth=min_total_depth,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SVIM:
        df_structural_variants = refine_svim_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            min_total_depth=min_total_depth,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.PBSV:
        df_structural_variants = refine_pbsv_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            keep_only_precise=keep_only_precise_sv,
            min_total_depth=min_total_depth,
            min_variant_reads_count=min_variant_reads_count
        )
    else:
        raise Exception(
            "Invalid value for 'variant_calling_method': %s. "
            "Allowed 'variant_calling_method' values are %s"
            % (variant_calling_method,
               ', '.join(f"'{item}'" for item in VariantCallingMethods.StructuralVariantCallingMethods.ALL))
        )

    # Filter out variants near the gapped regions.
    if df_gapped_regions is not None:
        df_structural_variants = remove_structural_variants_near_gapped_regions(
            df_structural_variants=df_structural_variants,
            df_gapped_regions=df_gapped_regions,
            gapped_regions_padding=gapped_regions_padding
        )
        logger.info(
            '%i variants after filtering out variants near the gapped regions.'
            % len(df_structural_variants)
        )

    # Filter out excluded variants.
    if df_structural_variants_to_exclude is not None:
        df_structural_variants = remove_structural_variants(
            df_structural_variants=df_structural_variants,
            df_structural_variants_to_exclude=df_structural_variants_to_exclude,
            exclude_variants_padding=exclude_variants_padding
        )
        logger.info('%i variants after filtering out excluded variants.' % len(df_structural_variants))

    logger.info('%i variants after refinement.' % len(df_structural_variants))
    return df_structural_variants


def run_exacto_refine_genomic_small_variants(
        df_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        variant_calling_method: str,
        is_tumor_normal_paired: bool,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding: int = GENOME_GAPPED_REGIONS_PADDING) -> pd.DataFrame:
    """
    Refines a set of variants and returns the refined set.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'chrom'
                                'pos'
                                'tumor_total_coverage'
                                'tumor_variant_reads_count'
    df_gapped_regions       :   DataFrame of gapped regions in the genome.
    variant_calling_method  :   Variant calling method.
    is_tumor_normal_paired  :   If true, variants were called using
                                tumor and matched normal pair.
    keep_only_chromosomes   :   List of chromosomes to keep.
                                Chromosomes not specified in this list
                                will be removed.
    keep_only_filter_values :   List of 'filter' values to keep.
                                'filter' values not specified in this list
                                will be removed.
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum variant reads count.
    gapped_regions_padding  :   Number of bases to pad gapped regions'
                                start and end positions.

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2:
        if is_tumor_normal_paired:
            df_variants = refine_gatk4_mutect2_tumor_normal_callset(
                df_variants=df_variants,
                keep_only_chromosomes=keep_only_chromosomes,
                keep_only_filter_values=keep_only_filter_values,
                min_total_depth=min_total_depth,
                min_variant_reads_count=min_variant_reads_count,
                min_normal_total_depth=min_total_depth
            )
        else:
            df_variants = refine_gatk4_mutect2_tumor_only_callset(
                df_variants=df_variants,
                keep_only_chromosomes=keep_only_chromosomes,
                keep_only_filter_values=keep_only_filter_values,
                min_total_depth=min_total_depth,
                min_variant_reads_count=min_variant_reads_count
            )
    elif variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT:
        df_variants = refine_deepvariant_callset(
            df_variants=df_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            min_total_depth=min_total_depth,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.STRELKA2:
        if is_tumor_normal_paired:
            df_variants = refine_strelka2_tumor_normal_callset(
                df_variants=df_variants,
                keep_only_chromosomes=keep_only_chromosomes,
                keep_only_filter_values=keep_only_filter_values,
                min_total_depth=min_total_depth,
                min_variant_reads_count=min_variant_reads_count,
                min_normal_total_depth=min_total_depth
            )
        else:
            df_variants = refine_strelka2_tumor_only_callset(
                df_variants=df_variants,
                keep_only_chromosomes=keep_only_chromosomes,
                keep_only_filter_values=keep_only_filter_values,
                min_total_depth=min_total_depth,
                min_variant_reads_count=min_variant_reads_count
            )
    else:
        raise Exception(
            "Invalid value for 'variant_calling_method': %s. "
            "Allowed 'variant_calling_method' values are %s"
            % (variant_calling_method,
               ', '.join(f"'{item}'" for item in VariantCallingMethods.SmallVariantCallingMethods.ALL))
        )

    # Filter out variants near the gapped regions.
    if df_gapped_regions is not None:
        df_variants = remove_small_variants_near_gapped_regions(
            df_variants=df_variants,
            df_gapped_regions=df_gapped_regions,
            gapped_regions_padding=gapped_regions_padding
        )
        logger.info(
            '%i variants after filtering out variants near the gapped regions.'
            % len(df_variants)
        )

    logger.info('%i variants after refinement.' % len(df_variants))
    return df_variants


def run_exacto_annotate_genomic_small_variants(
        df_small_variants: pd.DataFrame,
        annotation_source: str,
        df_gencode_genes: pd.DataFrame,
        df_gencode_exons: pd.DataFrame,
        ensembl_release: int = -1,
        perl_path: str = '',
        annovar_path: str = '',
        annovar_humandb_path: str = '',
        annovar_protocol: str = '',
        annovar_operation: str = '',
        annovar_genome_assembly: str = '',
        annovar_avinput_file: str = '',
        annovar_output_file: str = '') -> pd.DataFrame:
    """
    Annotates a set of variants and returns the annotated set.

    Parameters
    ----------
    df_small_variants       :   DataFrame of variants. Expected columns:
                                'chrom'
                                'pos'
    annotation_source       :   Annotation source ('ensembl' or 'gencode').
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Specify this if 'annotation_source' is 'gencode'.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Specify this if 'annotation_source' is 'gencode'.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'
    ensembl_release         :   Ensembl release version.
                                Specify this if 'annotation_source' is 'ensembl'.
    perl_path               :   perl path.
    annovar_path            :   ANNOVAR path.
    annovar_humandb_path    :   ANNOVAR humandb/ directory path.
    annovar_protocol        :   ANNOVAR protocol.
    annovar_operation       :   ANNOVAR operation.
    annovar_genome_assembly :   ANNOVAR genome assembly.
    annovar_avinput_file    :   ANNOVAR AVINPUT file.
    annovar_output_file     :   ANNOVAR output file.

    Returns
    -------
    DataFrame of annotated genomic structural variants.
    """
    if annotation_source == AnnotationSources.ENSEMBL:
        df_small_variants = annotate_small_variants_using_pyensembl(
            df_small_variants=df_small_variants,
            ensembl_release=ensembl_release
        )
    elif annotation_source == AnnotationSources.GENCODE:
        df_small_variants = annotate_small_variants_using_gencode(
            df_small_variants=df_small_variants,
            df_gencode_genes=df_gencode_genes,
            df_gencode_exons=df_gencode_exons
        )
    elif annotation_source == AnnotationSources.ANNOVAR:
        df_small_variants = annotate_small_variants_using_annovar(
            perl_path=perl_path,
            annovar_path=annovar_path,
            humandb_path=annovar_humandb_path,
            avinput_file=annovar_avinput_file,
            genome_assembly=annovar_genome_assembly,
            protocol=annovar_protocol,
            operation=annovar_operation,
            output_file=annovar_output_file
        )
    else:
        raise Exception(
            "Invalid value for 'annotation_source': %s. Allowed 'annotation_source' values are %s "
            % (annotation_source, ', '.join(AnnotationSources.ALL)))
    return df_small_variants


def run_exacto_annotate_genomic_structural_variants(
        df_structural_variants: pd.DataFrame,
        annotation_source: str,
        df_gencode_genes: pd.DataFrame,
        df_gencode_exons: pd.DataFrame,
        ensembl_release: int = -1) -> pd.DataFrame:
    """
    Annotates a set of variants and returns the annotated set.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of variants.
                                Expected columns:
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
                                'sv_type' (DEL, INS, INV, DUP, BND or TRA)
    annotation_source       :   Annotation source ('ensembl' or 'gencode').
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Specify this if 'annotation_source' is 'gencode'.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Specify this if 'annotation_source' is 'gencode'.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'
    ensembl_release         :   Ensembl release version.
                                Specify this if 'annotation_source' is 'ensembl'.

    Returns
    -------
    DataFrame of annotated genomic structural variants.
    """
    if len(df_structural_variants) == 0:
        logger.warning('DataFrame is empty. Returning without annotating.')
        return df_structural_variants
    if annotation_source == AnnotationSources.ENSEMBL:
        df_structural_variants = annotate_structural_variants_using_pyensembl(
            df_structural_variants=df_structural_variants,
            ensembl_release=ensembl_release
        )
    elif annotation_source == AnnotationSources.GENCODE:
        df_structural_variants = annotate_structural_variants_using_gencode(
            df_structural_variants=df_structural_variants,
            df_gencode_genes=df_gencode_genes,
            df_gencode_exons=df_gencode_exons
        )
    else:
        raise Exception(
            "Invalid value for 'annotation_source': %s. Allowed 'annotation_source' values are %s "
            % (annotation_source, ', '.join(AnnotationSources.ALL)))
    return df_structural_variants


def run_exacto_merge_genomic_structural_variants(
        list_df: List[pd.DataFrame],
        enforce_variant_type_matching: bool = True,
        max_clustering_distance: int = MAX_SV_CLUSTER_DISTANCE) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merges a list of structural variant DataFrames.

    Parameters
    ----------
    list_df                         :   List of DataFrames.
                                        Expected columns in each DataFrame:
                                        'variant_id'
                                        'chr_1'
                                        'pos_1'
                                        'chr_2'
                                        'pos_2'
                                        'sv_type'
                                        'variant_calling_method'
                                        'sequencing_platform'
    enforce_variant_type_matching    :  If true, sv_type must match for two SVs
                                        to be merged into one.
    max_clustering_distance          :  Maximum SV clustering distance (default: 10).

    Returns
    -------
    df_merged                       :   DataFrame of all variants from all DataFrames merged.
    df_merged_deduped               :   DataFrame of all variants from all DataFrames merged and deduped.
    """
    df_merged, df_merged_deduped = merge_structural_variants(
        list_df=list_df,
        enforce_variant_type_matching=enforce_variant_type_matching,
        max_clustering_distance=max_clustering_distance
    )
    return df_merged, df_merged_deduped


def run_exacto_merge_annotations(list_df: List[pd.DataFrame]) -> pd.DataFrame:
    """
    Merges a list of DataFrames.

    Parameters
    ----------
    list_df                             :   List of DataFrames.
                                            Expected columns in each DataFrame:
                                            'variant_id'

    Returns
    -------
    df_merged               :   DataFrame of all annotations from all DataFrames merged.
    """
    df_merged = merge_annotations(list_df=list_df)
    return df_merged


def run_exacto_merge_genomic_small_variants(
        list_df: List[pd.DataFrame],
        enforce_variant_type_matching: bool = True,
        max_clustering_distance: int = MAX_SMALL_VARIANT_CLUSTER_DISTANCE) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merges a list of small variant DataFrames.

    Parameters
    ----------
    list_df                             :   List of DataFrames.
                                            Expected columns in each DataFrame:
                                            'variant_id'
                                            'chrom'
                                            'pos'
                                            'variant_type'
                                            'variant_calling_method'
                                            'sequencing_platform'
    enforce_variant_type_matching       :   If true, sv_type must match for two SVs
                                            to be merged into one.
    max_clustering_distance             :   Maximum small variant clustering distance (default: 1).

    Returns
    -------
    df_merged               :   DataFrame of all variants from all DataFrames merged.
    df_merged_deduped       :   DataFrame of all variants from all DataFrames merged and deduped.
    """
    df_merged, df_merged_deduped = merge_small_variants(
        list_df=list_df,
        enforce_variant_type_matching=enforce_variant_type_matching,
        max_clustering_distance=max_clustering_distance
    )
    return df_merged, df_merged_deduped



def run_exacto_identify_rna_variants(bam_file: pysam.AlignmentFile,
                                     num_cores: int) -> pd.DataFrame:
    """
    Identifies RNA variants in a BAM file.

    Parameters
    ----------
    bam_file    :   Pysam AlignmentFile object of a BAM file.
    num_cores   :   Number of cores to use.

    Returns
    -------
    """
    bam_filename = str(bam_file.filename.decode())
    variant_callset = exactors.identify_rna_variants(bam_filename, num_cores)
    df_variants = pd.DataFrame({
        'chromosome': variant_callset.chromosomes,
        'position': variant_callset.positions,
        'read_id': variant_callset.read_ids,
        'variant_type': variant_callset.variant_types,
        'reference_allele': variant_callset.reference_alleles,
        'alternate_allele': variant_callset.alternate_alleles,
        'sequence': variant_callset.sequences,
        'variant_size': variant_callset.variant_sizes
    })
    return df_variants


def run_exacto_simulate_rna_variants(
        genome_fasta: pysam.FastaFile,
        df_genes: pd.DataFrame,
        df_transcripts: pd.DataFrame,
        df_exons: pd.DataFrame,
        df_target_regions: pd.DataFrame,
        df_herv_regions: pd.DataFrame,
        num_snv: int = SIMULATE_RNA_VARIANTS_NUM_SNV,
        num_insertion: int = SIMULATE_RNA_VARIANTS_NUM_INSERTION,
        num_deletion: int = SIMULATE_RNA_VARIANTS_NUM_DELETION,
        num_fusion: int = SIMULATE_RNA_VARIANTS_NUM_FUSION,
        num_inversion: int = SIMULATE_RNA_VARIANTS_NUM_INVERSION,
        num_herv: int = SIMULATE_RNA_VARIANTS_NUM_HERV,
        insertion_size_mean: int = SIMULATE_RNA_VARIANTS_INSERTION_SIZE_MEAN,
        insertion_size_stdev: int = SIMULATE_RNA_VARIANTS_INSERTION_SIZE_STDEV,
        deletion_size_mean: int = SIMULATE_RNA_VARIANTS_DELETION_MEAN,
        deletion_size_stdev: int = SIMULATE_RNA_VARIANTS_DELETION_STDEV,
        herv_solo_ltr_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_SOLO_LTR,
        herv_truncated_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_TRUNCATED,
        herv_chimeric_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_CHIMERIC,
        herv_chimeric_max_neighboring_distance: int = SIMULATE_RNA_VARIANTS_HERV_CHIMERIC_MAX_NEIGHBORING_DISTANCE,
        herv_full_length_proportion: float = SIMULATE_RNA_VARIANTS_HERV_PROPORTION_FULL_LENGTH,
        infinite_sites_assumption: bool = SIMULATE_RNA_VARIANTS_INFINITE_SITES_ASSUMPTION) -> Tuple[pd.DataFrame, List]:
    """
    Simulates RNA variants.

    Parameters
    ----------
    genome_fasta                            :   pysam.FastaFile object of reference genome.
    df_transcripts                          :   DataFrame of transcripts.
    df_exons                                :   DataFrame of exons.
    df_target_regions                       :   DataFrame of regions to simulate RNA variants.
    df_herv_regions                         :   HERV regions.
    num_snv                                 :   Number of SNVs to simulate.
    num_insertion                           :   Number of insertions to simulate.
    num_deletion                            :   Number of deletions to simulate.
    num_fusion                              :   Number of fusions to simulate.
    num_inversion                           :   Number of inversions to simulate.
    num_herv                                :   Number of HERVs to simulate.
    insertion_size_mean                     :   Mean value of insertion size.
    insertion_size_stdev                    :   Standard deviation of insertion size.
    deletion_size_mean                      :   Mean value of deletion size.
    deletion_size_stdev                     :   Standard deviation of deletion size.
    herv_solo_ltr_proportion                :   Proportion of expressed HERVs that only have solo LTR sequences.
    herv_truncated_proportion               :   Proportion of HERVs that are truncated.
    herv_chimeric_proportion                :   Proportion of HERVs that are chimeric (concatenation of neighboring HERVs).
    herv_chimeric_max_neighboring_distance  :   Maximum distance for two HERVs to be considered for simulation of a
                                                chimeric HERV.
    herv_full_length_proportion             :   Proportion of HERVs that are full-lengths.
    infinite_sites_assumption               :   If true, the simulation enforces the infinite sites assumption.

    Returns
    -------
    df_rna_variants                         :   DataFrame of RNA variants
    variant_transcript_sequences            :   List of variant transcript sequences
    """
    df_rna_variants, variant_transcript_sequences = simulate_rna_variants(**locals())
    return df_rna_variants, variant_transcript_sequences


def run_exacto_simulate_reads(sequences: List[Sequence],
                              num_gigabases: float,
                              read_length_mean: float,
                              read_length_stdev: float,
                              base_quality_mean: float,
                              base_quality_stdev: float) -> List[Read]:
    """
    Simulates sequencing reads.

    Parameters
    ----------
    sequences           :   List of instances of the class Sequence.
    num_gigabases       :   Number of gigabases to sequence.
    read_length_mean    :   Mean value of read length.
    read_length_stdev   :   Standard deviation of read length.
    base_quality_mean   :   Mean value of base quality.
    base_quality_stdev  :   Standard deviation of base quality.

    Returns
    -------
    reads               :   List of instances of the class Read.
    """
    return simulate_single_end_reads(sequences=sequences,
                                     num_bases=num_gigabases * 10e9,
                                     read_length_mean=read_length_mean,
                                     read_length_stdev=read_length_stdev,
                                     base_quality_mean=base_quality_mean,
                                     base_quality_stdev=base_quality_stdev)
