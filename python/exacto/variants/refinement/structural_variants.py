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
The purpose of this python3 script is to implement functions used to
refine structural variants.
"""


import pandas as pd
import numpy as np
import multiprocessing as mp
from typing import List
from ...default_parameters import *
from ...logging import get_logger
from ...utilities.pandas import *


logger = get_logger(__name__)


def remove_structural_variants_near_gapped_regions_work(
        df_structural_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame
    ):
    keep = []
    for index, row in df_structural_variants.iterrows():
        conditions = ((df_gapped_regions['chrom'] == row['chr_1']) &
                      (df_gapped_regions['start'] <= row['pos_1']) &
                      (df_gapped_regions['end'] >= row['pos_1'])) | \
                     ((df_gapped_regions['chrom'] == row['chr_2']) &
                      (df_gapped_regions['start'] <= row['pos_2']) &
                      (df_gapped_regions['end'] >= row['pos_2']))
        df_matched = df_gapped_regions[conditions]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df_structural_variants = df_structural_variants.loc[keep, :]
    return df_structural_variants


def remove_structural_variants_near_gapped_regions(
        df_structural_variants: pd.DataFrame,
        df_gapped_regions: pd.DataFrame,
        gapped_regions_padding: int = GENOME_GAPPED_REGIONS_PADDING,
        num_processes: int = NUM_PROCESSES_REFINE
    ) -> pd.DataFrame:
    """
    Removes structural variants with breakpoints near gapped regions.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of structural variants.
                                Expected columns:
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
    df_gapped_regions       :   DataFrame of gapped regions.
                                Expected columns:
                                'chrom'
                                'chromStart'
                                'chromEnd'
    gapped_regions_padding  :   Gapped regions padding.
    num_processes           :   Number of processes.

    Returns
    -------
    DataFrame of structural variants.
    """
    df_gapped_regions['start'] = df_gapped_regions.apply(
        lambda row: int(row.chromStart - gapped_regions_padding), axis=1
    )
    df_gapped_regions['end'] = df_gapped_regions.apply(
        lambda row: int(row.chromEnd + gapped_regions_padding), axis=1
    )
    list_df = np.array_split(df_structural_variants, num_processes)
    pool = mp.Pool(processes=num_processes)
    async_results = [pool.apply_async(remove_structural_variants_near_gapped_regions_work, args=(df_curr, df_gapped_regions)) for df_curr in list_df]
    pool.close()
    pool.join()
    df_structural_variants_list = [ar.get() for ar in async_results]
    df_structural_variants = pd.concat(df_structural_variants_list)
    return df_structural_variants


def remove_structural_variants_work(
        df_structural_variants: pd.DataFrame,
        df_structural_variants_to_exclude: pd.DataFrame
    ):
    keep = []
    for index, row in df_structural_variants.iterrows():
        conditions = \
            ((df_structural_variants_to_exclude['chr_1'] == row['chr_1']) &
             (df_structural_variants_to_exclude['pos_1_start'] <= row['pos_1']) &
             (df_structural_variants_to_exclude['pos_1_end'] >= row['pos_1']) &
             (df_structural_variants_to_exclude['sv_type'] == row['sv_type'])) | \
            ((df_structural_variants_to_exclude['chr_2'] == row['chr_1']) &
             (df_structural_variants_to_exclude['pos_2_start'] <= row['pos_1']) &
             (df_structural_variants_to_exclude['pos_2_end'] >= row['pos_1']) &
             (df_structural_variants_to_exclude['sv_type'] == row['sv_type'])) | \
            ((df_structural_variants_to_exclude['chr_1'] == row['chr_2']) &
             (df_structural_variants_to_exclude['pos_1_start'] <= row['pos_2']) &
             (df_structural_variants_to_exclude['pos_1_end'] >= row['pos_2']) &
             (df_structural_variants_to_exclude['sv_type'] == row['sv_type'])) | \
            ((df_structural_variants_to_exclude['chr_2'] == row['chr_2']) &
             (df_structural_variants_to_exclude['pos_2_start'] <= row['pos_2']) &
             (df_structural_variants_to_exclude['pos_2_end'] >= row['pos_2']) &
             (df_structural_variants_to_exclude['sv_type'] == row['sv_type']))
        df_matched = df_structural_variants_to_exclude[conditions]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df_structural_variants = df_structural_variants.loc[keep, :]
    return df_structural_variants


def remove_structural_variants(
        df_structural_variants: pd.DataFrame,
        df_structural_variants_to_exclude: pd.DataFrame,
        exclude_variants_padding: int = EXCLUDE_SV_PADDING,
        num_processes: int = NUM_PROCESSES_REFINE
    ) -> pd.DataFrame:
    """
    Removes structural variants with breakpoints near a list of structural variants to exclude.

    Parameters
    ----------
    df_structural_variants              :   DataFrame of structural variants.
                                            Expected columns:
                                            'chr_1'
                                            'pos_1'
                                            'chr_2'
                                            'pos_2'
                                            'sv_type'
    df_structural_variants_to_exclude   :   DataFrame of structural variants to exclude.
                                            Expected columns:
                                            'chr_1'
                                            'pos_1'
                                            'chr_2'
                                            'pos_2'
                                            'sv_type'
    exclude_variants_padding            :   Number of bases to pad breakpoints of
                                            structural variants in df_structural_variants_to_exclude
    num_processes                       :   Number of processes.

    Returns
    -------
    DataFrame of structural variants.
    """
    df_structural_variants_to_exclude['pos_1_start'] = df_structural_variants_to_exclude.apply(
        lambda row: int(row.pos_1 - exclude_variants_padding), axis=1
    )
    df_structural_variants_to_exclude['pos_1_end'] = df_structural_variants_to_exclude.apply(
        lambda row: int(row.pos_1 + exclude_variants_padding), axis=1
    )
    df_structural_variants_to_exclude['pos_2_start'] = df_structural_variants_to_exclude.apply(
        lambda row: int(row.pos_2 - exclude_variants_padding), axis=1
    )
    df_structural_variants_to_exclude['pos_2_end'] = df_structural_variants_to_exclude.apply(
        lambda row: int(row.pos_2 + exclude_variants_padding), axis=1
    )
    list_df = np.array_split(df_structural_variants, num_processes)
    pool = mp.Pool(processes=num_processes)
    async_results = [pool.apply_async(remove_structural_variants_work, args=(df_curr, df_structural_variants_to_exclude)) for df_curr in list_df]
    pool.close()
    pool.join()
    df_structural_variants_list = [ar.get() for ar in async_results]
    df_structural_variants = pd.concat(df_structural_variants_list)
    return df_structural_variants


def refine_sniffles2_sv_callset(
        df_structural_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        keep_only_precise: bool = KEEP_ONLY_PRECISE_SV,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines structural variants called using Sniffles2 and returns a DataFrame
    of refined structural variants.

    Parameters
    ----------
    df_structural_variants          :   DataFrame of structural variants.
                                        Expected columns:
                                        'chr_1'
                                        'chr_2'
                                        'filter'
                                        'is_precise'
                                        'total_depth'
                                        'variant_reads_count'
    keep_only_chromosomes           :   List of chromosomes to keep.
                                        Chromosomes not specified in this list
                                        will be filtered out.
    keep_only_filter_values         :   List of FILTER values to include.
    keep_only_precise               :   If true, only structural variants with
                                        'precise' breakpoints are kept.
    min_total_depth                 :   Minimum total depth (default: 7).
    min_variant_reads_count         :   Minimum number of variants (support) reads (default: 3).

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if len(keep_only_filter_values) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['filter'].isin(keep_only_filter_values)
        ]
    if keep_only_precise:
        df_structural_variants = df_structural_variants[
            df_structural_variants['is_precise']
        ]
    if len(keep_only_chromosomes) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['chr_1'].isin(keep_only_chromosomes) &
            df_structural_variants['chr_2'].isin(keep_only_chromosomes)
        ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'] >= min_total_depth
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'] >= min_variant_reads_count
    ]
    return df_structural_variants


def refine_cutesv_sv_callset(
        df_structural_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        keep_only_precise: bool = KEEP_ONLY_PRECISE_SV,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines structural variants called using cuteSV and returns a DataFrame
    of refined structural variants.

    Parameters
    ----------
    df_structural_variants          :   DataFrame of structural variants.
                                        Expected columns:
                                        'chr_1'
                                        'chr_2'
                                        'filter'
                                        'is_precise'
                                        'total_coverage'
                                        'variant_reads_count'
    keep_only_chromosomes           :   List of chromosomes to keep.
                                        Chromosomes not specified in this list
                                        will be filtered out.
    keep_only_filter_values         :   List of FILTER values to include.
    keep_only_precise               :   If true, only structural variants with
                                        'precise' breakpoints are kept.
    min_total_depth                 :   Minimum total depth (default: 7).
    min_variant_reads_count         :   Minimum number of variants (support) reads (default: 3).

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if len(keep_only_filter_values) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['filter'].isin(keep_only_filter_values)
        ]
    if keep_only_precise:
        df_structural_variants = df_structural_variants[
            df_structural_variants['is_precise'] == True
        ]
    if len(keep_only_chromosomes) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['chr_1'].isin(keep_only_chromosomes) &
            df_structural_variants['chr_2'].isin(keep_only_chromosomes)
        ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'] >= min_total_depth
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'] >= min_variant_reads_count
    ]
    return df_structural_variants


def refine_svim_sv_callset(
        df_structural_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines structural variants called using SVIM and returns a DataFrame
    of refined structural variants.

    Parameters
    ----------
    df_structural_variants          :   DataFrame of structural variants.
                                        Expected columns:
                                        'chr_1'
                                        'chr_2'
                                        'filter'
                                        'total_coverage'
                                        'variant_reads_count'
    keep_only_chromosomes           :   List of chromosomes to keep.
                                        Chromosomes not specified in this list
                                        will be filtered out.
    keep_only_filter_values         :   List of FILTER values to include.
    min_total_depth                 :   Minimum total depth (default: 7).
    min_variant_reads_count         :   Minimum number of variants (support) reads (default: 3).

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if len(keep_only_filter_values) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['filter'].isin(keep_only_filter_values)
        ]
    if len(keep_only_chromosomes) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['chr_1'].isin(keep_only_chromosomes) &
            df_structural_variants['chr_2'].isin(keep_only_chromosomes)
        ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'] >= min_total_depth
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'] >= min_variant_reads_count
    ]
    return df_structural_variants


def refine_pbsv_sv_callset(
        df_structural_variants: pd.DataFrame,
        keep_only_chromosomes: List[str] = [],
        keep_only_filter_values: List[str] = KEEP_ONLY_FILTER_VALUES,
        keep_only_precise: bool = KEEP_ONLY_PRECISE_SV,
        min_total_depth: int = MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count: int = MIN_GENOMIC_VARIANT_READS_COUNT
    ) -> pd.DataFrame:
    """
    Refines structural variants called using PBSV and returns a DataFrame
    of refined structural variants.

    Parameters
    ----------
    df_structural_variants          :   DataFrame of structural variants.
                                        Expected columns:
                                        'chr_1'
                                        'chr_2'
                                        'filter'
                                        'is_precise'
                                        'total_depth'
                                        'variant_reads_count'
    keep_only_chromosomes           :   List of chromosomes to keep.
                                        Chromosomes not specified in this list
                                        will be filtered out.
    keep_only_filter_values         :   List of FILTER values to include.
    keep_only_precise               :   If true, only structural variants with
                                        'precise' breakpoints are kept.
    min_total_depth                 :   Minimum total depth (default: 7).
    min_variant_reads_count         :   Minimum number of variants (support) reads (default: 3).

    Returns
    -------
    DataFrame of refined structural variants.
    """
    if len(keep_only_filter_values) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['filter'].isin(keep_only_filter_values)
        ]
    if keep_only_precise:
        df_structural_variants = df_structural_variants[
            df_structural_variants['is_precise'] == True
        ]
    if len(keep_only_chromosomes) > 0:
        df_structural_variants = df_structural_variants[
            df_structural_variants['chr_1'].isin(keep_only_chromosomes) &
            df_structural_variants['chr_2'].isin(keep_only_chromosomes)
        ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['total_depth'] >= min_total_depth
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'].map(is_safe_integer)
    ]
    df_structural_variants = df_structural_variants[
        df_structural_variants['variant_reads_count'] >= min_variant_reads_count
    ]
    return df_structural_variants
