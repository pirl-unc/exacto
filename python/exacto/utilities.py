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
The purpose of this python3 script is to implement general-purpose utility functions.
"""


import pandas as pd
from typing import Dict
from .constants import VariantCallingMethods
from .logging import get_logger


logger = get_logger(__name__)


def retrieve_with_default(dict, key, default_value, type):
    """
    Safely retrieves a value from a dictionary.

    Parameters
    ----------
    dict            :   Dictionary (or vector).
    key             :   Key.
    default_value   :   Default value.
    type            :   Type (str, int, float).

    Returns
    -------
    value           :   Value converted to 'type'.
                        If the conversion fails, the default value is returned.
    """
    # Try retrieving the value
    try:
        value = dict[key]
    except:
        value = default_value

    # Try converting the value to the desired type
    try:
        if type == str:
            value = str(value)
        elif type == int:
            value = int(value)
        elif type == float:
            value = float(value)
        else:
            value = default_value
    except:
        value = default_value
    return value


def get_typed_value(value, default_value, type):
    """
    Safely converts a value from a VCF row.

    Parameters
    ----------
    value           :   Value.
    default_value   :   Default value.
    type            :   Type (str, int, float).

    Returns
    -------
    value           :   Value converted to 'type'.
                        If the conversion fails, the default value is returned.
    """
    try:
        if type == str:
            value = str(value)
        elif type == int:
            value = int(value)
        elif type == float:
            value = float(value)
        elif type == bool:
            value = bool(value)
        else:
            value = default_value
    except:
        value = default_value
    return value


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


def get_variant_calling_method_attr_types(variant_calling_method: str) -> Dict:
    """
    Returns the attribute type dictionary for a given variant calling method.

    Parameters
    ----------
    variant_calling_method  :   Variant calling method.

    Returns
    -------
    attr_type_dict          :   Attribute-type dictionary.
    """
    if variant_calling_method == VariantCallingMethods.CUTESV:
        return VariantCallingMethods.AttributeTypes.CUTESV
    if variant_calling_method == VariantCallingMethods.DEEPVARIANT:
        return VariantCallingMethods.AttributeTypes.DEEPVARIANT
    if variant_calling_method == VariantCallingMethods.GATK4_MUTECT2:
        return VariantCallingMethods.AttributeTypes.GATK4_MUTECT2
    if variant_calling_method == VariantCallingMethods.PBSV:
        return VariantCallingMethods.AttributeTypes.PBSV
    if variant_calling_method == VariantCallingMethods.SNIFFLES2:
        return VariantCallingMethods.AttributeTypes.SNIFFLES2
    if variant_calling_method == VariantCallingMethods.STRELKA2:
        return VariantCallingMethods.AttributeTypes.STRELKA2
    if variant_calling_method == VariantCallingMethods.SVIM:
        return VariantCallingMethods.AttributeTypes.SVIM
