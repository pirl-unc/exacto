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
The purpose of this python3 script is to implement functions
related to merging structural variants.
"""


import pandas as pd
from typing import Tuple, List
from ...constants import *
from ...default_parameters import *
from ...logging import get_logger


logger = get_logger(__name__)


def merge_structural_variants(
        list_df: List[pd.DataFrame] = [],
        enforce_variant_type_matching: bool = True,
        max_clustering_distance: int = MAX_SV_CLUSTER_DISTANCE
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merges DataFrames of structural variants into one.

    Parameters
    ----------
    list_df                         :   List of DataFrames.
                                        Expected columns in each DataFrame:
                                        'chr_1'
                                        'pos_1'
                                        'chr_2'
                                        'pos_2'
                                        'sv_type'
                                        'variant_calling_method'
                                        'sequencing_platform'
    enforce_variant_type_matching   :   If true, sv_type must match for two SVs
                                        to be merged into one.
    max_clustering_distance         :   Maximum distance between methods to cluster (default: 10).

    Returns
    -------
    df_merged                       :   DataFrame of all variants from all DataFrames merged.
    df_merged_deduped               :   DataFrame of all variants from all DataFrames merged and deduped.
                                        The following columns are appended to each row:
                                        'matched_variant_ids'
                                        'variant_calling_methods'
                                        'variant_calling_methods_count'
    """
    df_merged = pd.DataFrame()
    logger.info('Started reading the DataFrames')
    i = 0
    shared_columns = []
    for df_curr in list_df:
        if len(df_curr) == 0:
            logger.warning("DataFrame is empty. Skipping merging this particular DataFrame.")
            continue
        if len(shared_columns) == 0:
            for curr_col in df_curr.columns.values.tolist():
                shared_columns.append(curr_col)
        else:
            shared_columns_temp = []
            for curr_col in df_curr.columns.values.tolist():
                if curr_col in shared_columns:
                    shared_columns_temp.append(curr_col)
            shared_columns = []
            for curr_col in shared_columns_temp:
                shared_columns.append(curr_col)
        df_merged = pd.concat([df_merged, df_curr])
        i += 1
    logger.info("Shared columns across all DataFrames:")
    logger.info(", ".join(shared_columns))
    df_merged['record_id'] = range(0, len(df_merged))
    logger.info('Finished reading the DataFrames')

    # Merge SV calls
    df_merged_deduped = pd.DataFrame()
    recorded_ids = set()
    n = len(df_merged)
    logger.info("%i SV calls to iterate" % n)
    for index, row in df_merged.iterrows():
        if row['record_id'] in recorded_ids:
            continue

        curr_variant_id = str(row['variant_id'])
        curr_chr_1 = str(row['chr_1'])
        curr_pos_1 = int(row['pos_1'])
        curr_chr_2 = str(row['chr_2'])
        curr_pos_2 = int(row['pos_2'])
        curr_sv_type = str(row['sv_type'])

        if enforce_variant_type_matching:
            if curr_sv_type == StructuralVariantTypes.INSERTION:
                query_sv_types = [StructuralVariantTypes.INSERTION,
                                  StructuralVariantTypes.DUPLICATION]
            elif curr_sv_type == StructuralVariantTypes.DELETION:
                query_sv_types = [StructuralVariantTypes.DELETION]
            elif curr_sv_type == StructuralVariantTypes.INVERSION:
                query_sv_types = [StructuralVariantTypes.INVERSION,
                                  StructuralVariantTypes.BREAKPOINT]
            elif curr_sv_type == StructuralVariantTypes.TRANSLOCATION:
                query_sv_types = [StructuralVariantTypes.TRANSLOCATION,
                                  StructuralVariantTypes.BREAKPOINT,
                                  StructuralVariantTypes.INVERSION]
            elif curr_sv_type == StructuralVariantTypes.BREAKPOINT:
                query_sv_types = [StructuralVariantTypes.BREAKPOINT,
                                  StructuralVariantTypes.INVERSION,
                                  StructuralVariantTypes.TRANSLOCATION]
            elif curr_sv_type == StructuralVariantTypes.DUPLICATION:
                query_sv_types = [StructuralVariantTypes.DUPLICATION,
                                  StructuralVariantTypes.INSERTION]
            else:
                raise Exception(
                    "Unknown SV type: %s. "
                    "Allowed SV types are: %s."
                    % (curr_sv_type, ', '.join(StructuralVariantTypes))
                )
            conditions = (
                (df_merged.sv_type.isin(query_sv_types)) &
                ((df_merged.chr_1 == curr_chr_1) & (df_merged.chr_2 == curr_chr_2) &
                 (df_merged.pos_1 >= curr_pos_1 - max_clustering_distance) &
                 (df_merged.pos_1 <= curr_pos_1 + max_clustering_distance) &
                 (df_merged.pos_2 >= curr_pos_2 - max_clustering_distance) &
                 (df_merged.pos_2 <= curr_pos_2 + max_clustering_distance)) |
                ((df_merged.chr_1 == curr_chr_2) & (df_merged.chr_2 == curr_chr_1) &
                 (df_merged.pos_1 >= curr_pos_2 - max_clustering_distance) &
                 (df_merged.pos_1 <= curr_pos_2 + max_clustering_distance) &
                 (df_merged.pos_2 >= curr_pos_1 - max_clustering_distance) &
                 (df_merged.pos_2 <= curr_pos_1 + max_clustering_distance))
            )
            df_matched = df_merged.loc[conditions,:]
        else:
            conditions = (
                ((df_merged.chr_1 == curr_chr_1) & (df_merged.chr_2 == curr_chr_2) &
                 (df_merged.pos_1 >= curr_pos_1 - max_clustering_distance) &
                 (df_merged.pos_1 <= curr_pos_1 + max_clustering_distance) &
                 (df_merged.pos_2 >= curr_pos_2 - max_clustering_distance) &
                 (df_merged.pos_2 <= curr_pos_2 + max_clustering_distance)) |
                ((df_merged.chr_1 == curr_chr_2) & (df_merged.chr_2 == curr_chr_1) &
                 (df_merged.pos_1 >= curr_pos_2 - max_clustering_distance) &
                 (df_merged.pos_1 <= curr_pos_2 + max_clustering_distance) &
                 (df_merged.pos_2 >= curr_pos_1 - max_clustering_distance) &
                 (df_merged.pos_2 <= curr_pos_1 + max_clustering_distance))
            )
            df_matched = df_merged.loc[conditions,:]

        for curr_record_id in df_matched.record_id.values.tolist():
            recorded_ids.add(curr_record_id)

        matched_variant_ids = []
        for curr_id in df_matched.variant_id.values.tolist():
            matched_variant_ids.append(curr_id)

        # Record data
        variant_calling_methods = df_matched.variant_calling_method.unique()
        variant_calling_methods_count = len(variant_calling_methods)
        df_curr = row.to_frame().T
        df_curr['matched_variant_ids'] = ','.join(matched_variant_ids)
        df_curr['variant_calling_methods'] = ','.join(variant_calling_methods)
        df_curr['variant_calling_methods_count'] = variant_calling_methods_count
        df_merged_deduped = pd.concat([df_merged_deduped, df_curr])

    logger.info("%i deduped variants in total" % len(df_merged_deduped))
    return df_merged, df_merged_deduped
