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
related to handling DeepVariant VCF files.
"""


import gzip
import pandas as pd
from collections import defaultdict
from .common import *
from ...logging import get_logger
from ...constants import *
from ...default_parameters import *


logger = get_logger(__name__)


def convert_deepvariant_vcf_to_dataframe(
        vcf_file: str,
        sequencing_platform: str,
        sample_id: str
    ) -> pd.DataFrame:
    """
    Convert a DeepVariant VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file            :   VCF file.
    sequencing_platform :   Sequencing platform.
    sample_id           :   Sample ID.

    Returns
    -------
    df                  :   DataFrame with the keys of
                            default_parameters.SMALL_VARIANT_ATTRIBUTES
                            as the columns
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['sample_id'] = sample_id
        curr_row['variant_id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['quality_score'] = safely_retrieve_value(dict=row, key='QUAL', default_value=curr_row['quality_score'], type='float')

        # Make sure 'chr' is in chrom
        if 'chr' not in curr_row['chrom'] and 'CHR' not in curr_row['chrom']:
            curr_row['chrom'] = 'chr' + curr_row['chrom']

        if len(curr_row['ref']) == 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SINGLE_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = curr_row['alt'].upper()
            curr_row['variant_size'] = 1
        elif len(curr_row['ref']) == 1 and len(curr_row['alt']) > 1:
            if ',' in curr_row['alt']:
                curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
                curr_row['variant_sequence'] = curr_row['alt'].upper()
            else:
                curr_row['variant_type'] = SmallVariantTypes.SMALL_INSERTION
                curr_row['variant_sequence'] = curr_row['alt'][1:].upper()
                curr_row['variant_size'] = len(curr_row['alt'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SMALL_DELETION
            curr_row['variant_sequence'] = curr_row['ref'][1:].upper()
            curr_row['variant_size'] = len(curr_row['ref'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) > 1:
            curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = curr_row['alt'].upper()
            curr_row['variant_size'] = len(curr_row['alt'])
        else:
            logger.warning('Unknown variant type. REF: %s. ALT: %s' %
                           (curr_row['ref'], curr_row['alt']))

        # Extract FORMAT
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = safely_retrieve_value(dict=sample, key=format.index('GT'), default_value=curr_row['genotype'], type='str')
        if 'GQ' in format:
            curr_row['genotype_quality'] = safely_retrieve_value(dict=sample, key=format.index('GQ'), default_value=curr_row['genotype_quality'], type='str')
        if 'DP' in format:
            curr_row['total_depth'] = safely_retrieve_value(dict=sample, key=format.index('DP'), default_value=curr_row['total_depth'], type='int')
        if 'AD' in format:
            curr_ad = safely_retrieve_value(dict=sample, key=format.index('AD'), default_value='unknown', type='str').split(',')
            try:
                curr_row['reference_reads_count'] = int(curr_ad[0])
            except:
                pass
            try:
                if curr_row['variant_type'] == SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT:
                    curr_row['variant_reads_count'] = ','.join(curr_ad[1:])
                else:
                    curr_row['variant_reads_count'] = int(curr_ad[1])
            except:
                pass
        if 'VAF' in format:
            if curr_row['variant_type'] == SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT:
                curr_row['variant_allele_fraction'] = safely_retrieve_value(dict=sample, key=format.index('VAF'), default_value=curr_row['variant_allele_fraction'], type='str')
            else:
                curr_row['variant_allele_fraction'] = safely_retrieve_value(dict=sample, key=format.index('VAF'), default_value=curr_row['variant_allele_fraction'], type='float')
        if 'PL' in format:
            curr_row['phred_scale_genotype_likelihoods'] = safely_retrieve_value(dict=sample, key=format.index('PL'), default_value=curr_row['phred_scale_genotype_likelihoods'], type='str')

        # Update total_depth if it is currently unknown but can be inferred
        if type(curr_row['variant_reads_count']) == int and \
            type(curr_row['reference_reads_count']) == int and \
            curr_row['total_depth'] == SMALL_VARIANT_ATTRIBUTES['total_depth']:
            curr_row['total_depth'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']

        # Update variant_allele_fraction if it is currently unknown
        if type(curr_row['variant_reads_count']) == int and \
            type(curr_row['total_depth']) == int and \
            curr_row['variant_allele_fraction'] == SMALL_VARIANT_ATTRIBUTES['variant_allele_fraction']:
            if curr_row['total_depth'] > 0:
                curr_row['variant_allele_fraction'] = float(curr_row['variant_reads_count']) / float(curr_row['total_depth'])

        # Update ID
        if curr_row['variant_id'] == '.':
            curr_row['variant_id'] = VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT + '.' + \
                             curr_row['variant_type'] + '.' + \
                             str(curr_idx)
            curr_idx += 1

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df
