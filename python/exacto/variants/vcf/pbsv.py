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
related to handling pbsv VCF files.
"""


import gzip
import pandas as pd
from collections import defaultdict
from .common import *
from ...logging import get_logger
from ...constants import *
from ...default_parameters import *


logger = get_logger(__name__)


def convert_pbsv_vcf_to_dataframe(
        vcf_file: str,
        sequencing_platform: str,
        sample_id: str
    ) -> pd.DataFrame:
    """
    Convert a PBSV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.
    sample_id               :   Sample ID.

    Returns
    -------
    df                      :   DataFrame with the keys of
                                default_parameters.STRUCTURAL_VARIANT_ATTRIBUTES
                                as the columns
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    curr_idx = 1
    included_mate_ids = set()
    for row in df_vcf.to_dict('records'):
        curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
        curr_row['sample_id'] = sample_id
        curr_row['variant_id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.StructuralVariantCallingMethods.PBSV
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chr_1'] = str(row['CHROM'])
        curr_row['chr_2'] = str(row['CHROM'])
        curr_row['pos_1'] = int(row['POS'])
        curr_row['ref'] = str(row['REF']).upper()
        curr_row['alt'] = str(row['ALT']).upper()
        curr_row['filter'] = str(row['FILTER'])
        curr_row['is_precise'] = True
        if row['QUAL'] != '.':
            curr_row['quality_score'] = float(row['QUAL'])

        # Extract INFO
        info = str(row['INFO']).split(';')
        mate_id = ''
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'PRECISE':
                curr_row['is_precise'] = True
            if curr_info_elements[0] == 'IMPRECISE':
                curr_row['is_precise'] = False
            if curr_info_elements[0] == 'END':
                curr_row['pos_2'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'SVTYPE':
                curr_row['sv_type'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'SVLEN':
                curr_row['sv_size'] = abs(int(curr_info_elements[1]))
            if curr_info_elements[0] == 'CIPOS':
                curr_row['ci_pos'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'MATEID':
                mate_id = str(curr_info_elements[1])
            if curr_info_elements[0] == 'SVANN':
                curr_row['repeat_annotation'] = str(curr_info_elements[1])

        if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT:
            # Check if current ID has been included
            if mate_id in included_mate_ids:
                continue
            included_mate_ids.add(curr_row['variant_id'])

        # Update chromosome 2 for 'BND'
        if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT:
            curr_id = curr_row['variant_id'].split("-")[1]
            curr_id = curr_id.split(":")[0]
            curr_row['chr_2'] = str(curr_id)

        # Update position 2 for 'BND'
        if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT:
            curr_id = curr_row['variant_id'].split("-")[1]
            curr_id = curr_id.split(":")[1]
            curr_row['pos_2'] = int(curr_id)

        # Make sure 'chr' is in chr_1 and chr_2
        if ('chr' not in curr_row['chr_1']) and ('CHR' not in curr_row['chr_2']) and (curr_row['chr_1'] != ''):
            curr_row['chr_1'] = 'chr' + curr_row['chr_1']
        if ('chr' not in curr_row['chr_2']) and ('CHR' not in curr_row['chr_2']) and (curr_row['chr_2'] != ''):
            curr_row['chr_2'] = 'chr' + curr_row['chr_2']

        # Update SV size for 'BND'
        if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT and curr_row['chr_1'] == curr_row['chr_2']:
            curr_row['sv_size'] = abs(curr_row['pos_2'] - curr_row['pos_1'])

        # Update strand for 'BND'
        if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT:
            alt_val = curr_row['alt']
            if (alt_val[0:2] == 'A[') or (alt_val[0:2] == 'C[') or (alt_val[0:2] == 'T[') or (alt_val[0:2] == 'G['):
                curr_row['strand'] = '+-'
            elif (alt_val[0:2] == 'A]') or (alt_val[0:2] == 'C]') or (alt_val[0:2] == 'T]') or (alt_val[0:2] == 'G]'):
                curr_row['strand'] = '-+'
            elif (alt_val[-2:] == ']A') or (alt_val[-2:] == ']C') or (alt_val[-2:] == ']T') or (alt_val[-2:] == ']G'):
                curr_row['strand'] = '+-'
            elif (alt_val[-2:] == '[A') or (alt_val[-2:] == '[C') or (alt_val[-2:] == '[T') or (alt_val[-2:] == '[G'):
                curr_row['strand'] = '-+'
            else:
                curr_row['strand'] = ''
        else:
            curr_row['strand'] = ''

        # Update insertion sequence
        if curr_row['sv_type'] == StructuralVariantTypes.INSERTION:
            curr_row['insertion_sequence'] = curr_row['alt']
        else:
            curr_row['insertion_sequence'] = ''

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = safely_retrieve_value(dict=sample, key=format.index('GT'), default_value=curr_row['genotype'], type='str')
        if 'DP' in format:
            curr_row['total_depth'] = safely_retrieve_value(dict=sample, key=format.index('DP'), default_value=curr_row['total_depth'], type='int')
        if 'AD' in format:
            try:
                curr_row['reference_reads_count'] = int(safely_retrieve_value(dict=sample, key=format.index('AD'), default_value=curr_row['reference_reads_count'], type='str').split(',')[0])
            except:
                pass
            try:
                curr_row['variant_reads_count'] = int(safely_retrieve_value(dict=sample, key=format.index('AD'), default_value=curr_row['variant_reads_count'], type='str').split(',')[1])
            except:
                pass
        if 'SAC' in format:
            curr_row['strand_reads'] = safely_retrieve_value(dict=sample, key=format.index('SAC'), default_value=curr_row['strand_reads'], type='str')

        # Update total_depth if it is currently unknown but can be inferred
        if type(curr_row['variant_reads_count']) == int and \
            type(curr_row['reference_reads_count']) == int and \
            curr_row['total_depth'] == STRUCTURAL_VARIANT_ATTRIBUTES['total_depth']:
            curr_row['total_depth'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']

        # Update variant_allele_fraction if it is currently unknown
        if type(curr_row['variant_reads_count']) == int and \
            type(curr_row['total_depth']) == int and \
            curr_row['variant_allele_fraction'] == STRUCTURAL_VARIANT_ATTRIBUTES['variant_allele_fraction']:
            if curr_row['total_depth'] > 0:
                curr_row['variant_allele_fraction'] = float(curr_row['variant_reads_count']) / float(curr_row['total_depth'])

        # Update ID
        if curr_row['variant_id'] == '.':
            curr_row['variant_id'] = VariantCallingMethods.StructuralVariantCallingMethods.PBSV + '.' + \
                             curr_row['variant_type'] + '.' + \
                             str(curr_idx)
            curr_idx += 1

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df
