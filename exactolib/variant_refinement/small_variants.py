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
The purpose of this python3 script is to implement functions that are used
to refine SNV and INDEL variants.
"""


import os
import pandas as pd
import subprocess as sp
from typing import Tuple, List
from ..logging import get_logger
from ..constants import *
from ..default_parameters import *
from ..utilities.vcf_utils import *
from ..utilities.pandas_utils import *


logger = get_logger(__name__)


def remove_small_variants_near_gapped_regions(
        df_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        gapped_regions_padding: int = GENOME_GAPPED_REGIONS_PADDING) -> pd.DataFrame:
    """
    Removes variants with breakpoints near gapped regions.

    Parameters
    ----------
    df_variants  :              DataFrame of structural variants.
                                Expected columns:
                                'chrom'
                                'pos'
    df_gapped_regions       :   DataFrame of gapped regions.
                                Expected columns:
                                'chrom'
                                'chromStart'
                                'chromEnd'
    gapped_regions_padding  :   Gapped regions padding.

    Returns
    -------
    DataFrame of variants.
    """
    df_gapped_regions['start'] = df_gapped_regions.apply(
        lambda row: int(row.chromStart - gapped_regions_padding), axis=1
    )
    df_gapped_regions['end'] = df_gapped_regions.apply(
        lambda row: int(row.chromEnd + gapped_regions_padding), axis=1
    )
    keep = []
    for index, row in df_variants.iterrows():
        # Check if variant position falls inside the (padded) gap region
        conditions = ((df_gapped_regions['chrom'] == row['chrom']) &
                      (df_gapped_regions['start'] <= row['pos']) &
                      (df_gapped_regions['end'] >= row['pos']))
        df_matched = df_gapped_regions[conditions]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df_variants = df_variants.loc[keep,:]
    return df_variants


def refine_gatk4_mutect2_tumor_normal_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_coverage: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT) -> pd.DataFrame:
    """
    Refines small variants called using GATK4-Mutect2 and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'tumor_total_coverage'
                                'tumor_variant_reads_count'
                                'normal_total_coverage'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_coverage      :   Minimum total coverage.
    min_variant_reads_count :   Minimum number of variants (support) reads.

    Returns
    -------
    DataFrame
    """
    if len(keep_only_filter_values) > 0:
        df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
    if len(keep_only_chromosomes) > 0:
        df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
    df_variants = df_variants[df_variants['tumor_total_coverage'] >= min_total_coverage]
    df_variants = df_variants[df_variants['tumor_variant_reads_count'] >= min_variant_reads_count]
    df_variants = df_variants[df_variants['normal_total_coverage'] >= min_total_coverage]
    return df_variants


def refine_deepvariant_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_coverage: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT) -> pd.DataFrame:
    """
    Refines small variants called using DeepVariant and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'variant_reads_count'
                                'total_read_depth'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_coverage      :   Minimum total coverage.
    min_variant_reads_count :   Minimum number of variants (support) reads.

    Returns
    -------
    DataFrame
    """
    if len(keep_only_filter_values) > 0:
        df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
    if len(keep_only_chromosomes) > 0:
        df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
    df_variants = df_variants[df_variants['total_coverage'].map(is_safe_integer)]
    df_variants = df_variants[df_variants['total_coverage'] >= min_total_coverage]
    df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
    df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
    return df_variants

