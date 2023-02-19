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
import numpy as np
import subprocess as sp
import multiprocessing as mp
from typing import Tuple, List
from ...logging import get_logger
from ...constants import *
from ...default_parameters import *
from ...utilities.pandas import *


logger = get_logger(__name__)


def remove_small_variants_near_gapped_regions_work(
        df_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame
    ):
    if len(df_variants) == 0:
        return pd.DataFrame()

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


def remove_small_variants_near_gapped_regions(
        df_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        padding: int = GENOME_GAPPED_REGIONS_PADDING,
        num_processes=1
    ) -> pd.DataFrame:
    """
    Removes variants with breakpoints near gapped regions.

    Parameters
    ----------
    df_variants         :   DataFrame of structural variants.
                            Expected columns:
                            'chrom'
                            'pos'
    df_gapped_regions   :   DataFrame of gapped regions.
                            Expected columns:
                            'chrom'
                            'chromStart'
                            'chromEnd'
    padding             :   Gapped regions padding.
    num_processes       :   Number of processes.

    Returns
    -------
    DataFrame of variants.
    """
    df_gapped_regions['start'] = df_gapped_regions.apply(
        lambda row: int(row.chromStart - padding), axis=1
    )
    df_gapped_regions['end'] = df_gapped_regions.apply(
        lambda row: int(row.chromEnd + padding), axis=1
    )
    list_df = np.array_split(df_variants, num_processes)
    pool = mp.Pool(processes=num_processes)
    async_results = [pool.apply_async(remove_small_variants_near_gapped_regions_work, args=(df_curr, df_gapped_regions)) for df_curr in list_df]
    pool.close()
    pool.join()
    df_variants_list = [ar.get() for ar in async_results]
    df_variants = pd.concat(df_variants_list)
    return df_variants


def remove_small_variants_work(
        df_variants: pd.DataFrame,
        df_exclude_snv_indel: pd.DataFrame,
        enforce_variant_type_check: bool
    ):
    if len(df_variants) == 0:
        return pd.DataFrame()

    keep = []
    for index, row in df_variants.iterrows():
        # Check if variant should be removed
        if enforce_variant_type_check:
            conditions = ((df_exclude_snv_indel['chrom'] == row['chrom']) &
                          (df_exclude_snv_indel['pos'] == row['pos']) &
                          (df_exclude_snv_indel['variant_type'] == row['variant_type']) &
                          (df_exclude_snv_indel['variant_sequence'] == row['variant_sequence']))
        else:
            conditions = ((df_exclude_snv_indel['chrom'] == row['chrom']) &
                          (df_exclude_snv_indel['pos'] == row['pos']))

        df_matched = df_exclude_snv_indel[conditions]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df_variants = df_variants.loc[keep, :]
    return df_variants


def remove_small_variants(
        df_variants: pd.DataFrame,
        df_exclude_snv_indel: pd.DataFrame,
        enforce_variant_type_check: bool = True,
        num_processes=1
    ) -> pd.DataFrame:
    """
    Exclude small variants that appear in a given list.

    Parameters
    ----------
    df_variants                 :   DataFrame of structural variants.
                                    Expected columns:
                                    'chrom'
                                    'pos'
                                    'variant_type'
                                    'variant_sequence'
    df_exclude_snv_indel        :   DataFrame of SNVs and INDELs to exclude.
                                    Expected columns:
                                    'chrom'
                                    'pos'
                                    'variant_type'
                                    'variant_sequence'
    enforce_variant_type_check  :   If True, checks for the following:
                                    'chrom', 'pos', 'variant_type', 'variant_sequence'
                                    If False, checks for the following:
                                    'chrom', 'pos'
    num_processes               :   Number of processes.

    Returns
    -------
    DataFrame of variants.
    """
    list_df = np.array_split(df_variants, num_processes)
    pool = mp.Pool(processes=num_processes)
    async_results = [pool.apply_async(remove_small_variants_work, args=(df_curr, df_exclude_snv_indel, enforce_variant_type_check)) for df_curr in list_df]
    pool.close()
    pool.join()
    df_variants_list = [ar.get() for ar in async_results]
    df_variants = pd.concat(df_variants_list)
    return df_variants


def refine_gatk4_mutect2_tumor_normal_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT,
        min_normal_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH
    ) -> pd.DataFrame:
    """
    Refines small variants called using GATK4-Mutect2 and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'total_depth'
                                'variant_reads_count'
                                'normal_total_depth'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum number of variants (support) reads.
    min_normal_total_depth  :   Minimum normal total depth.

    Returns
    -------
    DataFrame
    """
    try:
        if len(keep_only_filter_values) > 0:
            df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
        if len(keep_only_chromosomes) > 0:
            df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
        df_variants = df_variants[df_variants['total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['total_depth'] >= min_total_depth]
        df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
        df_variants = df_variants[df_variants['normal_total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['normal_total_depth'] >= min_normal_total_depth]
        return df_variants
    except:
        return pd.DataFrame(columns=SMALL_VARIANT_ATTRIBUTES.keys())


def refine_gatk4_mutect2_tumor_only_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines small variants called using GATK4-Mutect2 and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'total_depth'
                                'variant_reads_count'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum number of variants (support) reads.

    Returns
    -------
    DataFrame
    """
    try:
        if len(keep_only_filter_values) > 0:
            df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
        if len(keep_only_chromosomes) > 0:
            df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
        df_variants = df_variants[df_variants['total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['total_depth'] >= min_total_depth]
        df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
        return df_variants
    except:
        return pd.DataFrame(columns=SMALL_VARIANT_ATTRIBUTES.keys())


def refine_strelka2_tumor_normal_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT,
        min_normal_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH
    ) -> pd.DataFrame:
    """
    Refines small variants called using Strelka2 and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'total_depth'
                                'variant_reads_count'
                                'normal_total_depth'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum number of variants (support) reads.
    min_normal_total_depth  :   Minimum normal total depth.

    Returns
    -------
    DataFrame
    """
    try:
        if len(keep_only_filter_values) > 0:
            df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
        if len(keep_only_chromosomes) > 0:
            df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
        df_variants = df_variants[df_variants['total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['total_depth'] >= min_total_depth]
        df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
        df_variants = df_variants[df_variants['normal_total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['normal_total_depth'] >= min_normal_total_depth]
        return df_variants
    except:
        return pd.DataFrame(columns=SMALL_VARIANT_ATTRIBUTES.keys())


def refine_strelka2_tumor_only_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines small variants called using Strelka2 and returns a DataFrame
    of refined variants.

    Parameters
    ----------
    df_variants             :   DataFrame of variants.
                                Expected columns:
                                'filter'
                                'chrom'
                                'total_depth'
                                'variant_reads_count'
    keep_only_chromosomes   :   List of chromosomes to keep.
    keep_only_filter_values :   List of FILTER values to include.
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum number of variants (support) reads.

    Returns
    -------
    DataFrame
    """
    try:
        if len(keep_only_filter_values) > 0:
            df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
        if len(keep_only_chromosomes) > 0:
            df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
        df_variants = df_variants[df_variants['total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['total_depth'] >= min_total_depth]
        df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
        return df_variants
    except:
        return pd.DataFrame(columns=SMALL_VARIANT_ATTRIBUTES.keys())


def refine_deepvariant_callset(
        df_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
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
    min_total_depth         :   Minimum total depth.
    min_variant_reads_count :   Minimum number of variants (support) reads.

    Returns
    -------
    DataFrame
    """
    try:
        if len(keep_only_filter_values) > 0:
            df_variants = df_variants[df_variants['filter'].isin(keep_only_filter_values)]
        if len(keep_only_chromosomes) > 0:
            df_variants = df_variants[df_variants['chrom'].isin(keep_only_chromosomes)]
        df_variants = df_variants[df_variants['total_depth'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['total_depth'] >= min_total_depth]
        df_variants = df_variants[df_variants['variant_reads_count'].map(is_safe_integer)]
        df_variants = df_variants[df_variants['variant_reads_count'] >= min_variant_reads_count]
        return df_variants
    except:
        return pd.DataFrame(columns=SMALL_VARIANT_ATTRIBUTES.keys())
