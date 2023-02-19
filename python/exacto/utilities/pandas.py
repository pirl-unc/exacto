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
The purpose of this python3 script is to implement functions related to
parse pandas DataFrames.
"""


import pandas as pd


def is_safe_integer(x):
    try:
        int(x)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


def overlaps_any(
        df: pd.DataFrame,
        chrom: str,
        start: int,
        end: int
    ) -> bool:
    """
    Checks if a particular region overlaps with any regions in a DataFrame.

    Parameters
    ----------
    df          :   DataFrame with the following columns:
                    'chr_1', 'pos_1', 'chr_2', 'pos_2'
    chrom       :   Chromosome.
    start       :   Start position.
    end         :   End position.

    Returns
    -------
    True or False
    """
    # De Morgan's law on checking for non-overlapping regions
    df_matched = df.loc[
        (df['chr_1'] == chrom) &
        (df['chr_2'] == chrom) &
        (df['pos_2'] >= start) &
        (df['pos_1'] <= end),:
    ]
    if len(df_matched) > 0:
        return True
    else:
        return False
