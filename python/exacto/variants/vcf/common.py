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
The purpose of this python3 script is to implement common functions
related to handling VCF files.
"""


import gzip
import pandas as pd
from collections import defaultdict
from ...logging import get_logger
from ...constants import *
from ...default_parameters import *


logger = get_logger(__name__)


def safely_retrieve_value(
        dict,
        key,
        default_value,
        type
    ):
    """
    Safely converts a value from a VCF row.

    Parameters
    ----------
    dict            :   Dictionary (or vector).
    key             :   Key.
    default_value   :   Default value.
    type            :   Type ('str', 'int', 'float).

    Returns
    -------
    value           :   Value converted to 'type'.
                        If the conversion fails, the default value is returned.
    """
    # Step 1. Retrieve the value
    try:
        value = dict[key]
    except:
        value = default_value

    # Step 2. Convert value to the desired type
    try:
        if type == 'str':
            value = str(value)
        elif type == 'int':
            value = int(value)
        elif type == 'float':
            value = float(value)
        else:
            value = default_value
    except:
        value = default_value
    return value


def read_vcf_file(
        vcf_file: str
    ) -> pd.DataFrame:
    """
    Reads a VCF file and returns a DataFrame.

    Parameters
    ----------
    vcf_file    :   VCF file.

    Returns
    -------
    df_vcf      :   DataFrame of variants.
    """
    vcf_names = []
    is_gzipped = False
    if vcf_file.endswith(".gz"):
        is_gzipped = True
        with gzip.open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break
    else:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#CHROM"):
                    vcf_names = line.split('\t')
                    break

    vcf_names = [i.replace('\n', '') for i in vcf_names]
    vcf_names = ['CHROM' if i == '#CHROM' else i for i in vcf_names]

    if is_gzipped:
        df_vcf = pd.read_csv(vcf_file,
                             compression='gzip',
                             comment='#',
                             delim_whitespace=True,
                             header=None,
                             names=vcf_names)
    else:
        df_vcf = pd.read_csv(vcf_file,
                             comment='#',
                             delim_whitespace=True,
                             header=None,
                             names=vcf_names)
    return df_vcf
