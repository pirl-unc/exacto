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


import pysam
import pandas as pd
from typing import Tuple, List
from .constants import *
from .default_parameters import *
from .logging import get_logger
from .variant_refinement.structural_variants import *
from .variant_refinement.small_variants import *
from .variant_annotation.ensembl import *
from .variant_annotation.gencode import *
from .utilities.merging import *
from .utilities.gencode import *


logger = get_logger(__name__)


def run_exacto_refine_genomic_structural_variants(
        df_structural_variants: pd.DataFrame,
        df_structural_variants_to_exclude: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        variant_calling_method: str,
        keep_only_precise_sv: bool,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_coverage: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
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
    min_total_coverage                  :   Minimum total coverage.
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
            min_total_coverage=min_total_coverage,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.CUTESV:
        df_structural_variants = refine_cutesv_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            keep_only_precise=keep_only_precise_sv,
            min_total_coverage=min_total_coverage,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SVIM:
        df_structural_variants = refine_svim_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            min_total_coverage=min_total_coverage,
            min_variant_reads_count=min_variant_reads_count
        )
    elif variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.PBSV:
        df_structural_variants = refine_pbsv_sv_callset(
            df_structural_variants=df_structural_variants,
            keep_only_chromosomes=keep_only_chromosomes,
            keep_only_filter_values=keep_only_filter_values,
            keep_only_precise=keep_only_precise_sv,
            min_total_coverage=min_total_coverage,
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
        tumor_normal_paired: bool,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_coverage: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
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
    tumor_normal_paired     :   If true, variants were called using
                                tumor and matched normal pair.
    keep_only_chromosomes   :   List of chromosomes to keep.
                                Chromosomes not specified in this list
                                will be removed.
    keep_only_filter_values :   List of 'filter' values to keep.
                                'filter' values not specified in this list
                                will be removed.
    min_total_coverage      :   Minimum total coverage.
    min_variant_reads_count :   Minimum variant reads count.
    gapped_regions_padding  :   Number of bases to pad gapped regions'
                                start and end positions.

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2:
        if tumor_normal_paired:
            df_variants = refine_gatk4_mutect2_tumor_normal_callset(
                df_variants=df_variants,
                keep_only_chromosomes=keep_only_chromosomes,
                keep_only_filter_values=keep_only_filter_values,
                min_total_coverage=min_total_coverage,
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


def run_exacto_annotate_genomic_structural_variants(
        df_structural_variants: pd.DataFrame,
        annotation_source: str,
        df_gencode_genes: pd.DataFrame,
        df_gencode_exons: pd.DataFrame,
        ensembl_release: int) -> pd.DataFrame:
    """
    Annotates a set of variants and returns the annotated set.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of variants.
                                Expected columns:
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
        max_sv_cluster_distance: int = MAX_SV_CLUSTER_DISTANCE) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merges a list of structural variant DataFrames.

    Parameters
    ----------
    list_df                 :   List of DataFrames.
                                Expected columns in each DataFrame:
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
                                'sv_type'
                                'variant_calling_method'
                                'sequencing_platform'
    max_sv_cluster_distance :   Maximum SV clustering distance.

    Returns
    -------
    df_merged               :   DataFrame of all variants from all DataFrames merged.
    df_merged_deduped       :   DataFrame of all variants from all DataFrames merged and deduped.
    """
    df_merged, df_merged_deduped = merge_structural_variants(
        list_df=list_df,
        max_sv_cluster_distance=max_sv_cluster_distance
    )
    return df_merged, df_merged_deduped


def run_exacto_simulate_variants(nucleic_acid_type: str,
                                 fasta: pysam.FastaFile,
                                 num_snv: int = SIMULATE_NUM_SNV,
                                 num_insertion: int = SIMULATE_NUM_INSERTION,
                                 num_deletion: int = SIMULATE_NUM_DELETION) -> pd.DataFrame:
    df = pd.DataFrame()
    if nucleic_acid_type == NucleicAcidTypes.DNA:
        for chrom, size in zip(fasta.references, fasta.lengths):
            print(chrom, size)
    return df
