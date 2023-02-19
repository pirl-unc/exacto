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
related to merging annotations.
"""


import pandas as pd
from typing import Tuple, List
from ...constants import *
from ...default_parameters import *
from ...logging import get_logger


logger = get_logger(__name__)


def merge_annotations(
        list_df: List[pd.DataFrame] = []
    ) -> pd.DataFrame:
    """
    Merges annotations into one DataFrame.

    Parameters
    ----------
    list_df         :   List of DataFrames.
                        Expected columns in each DataFrame:
                        'variant_id'
    Returns
    -------
    df_merged       :   DataFrame of all annotations.
    """
    df_merged = list_df[0]
    for i in range(1, len(list_df)):
        df_merged = pd.concat([df_merged, list_df[i]], axis=1, join="inner")
    return df_merged
