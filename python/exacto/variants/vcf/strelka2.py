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
related to handling Strelka2 VCF files.
"""


import gzip
import pandas as pd
from collections import defaultdict
from .common import *
from ...logging import get_logger
from ...constants import *
from ...default_parameters import *


logger = get_logger(__name__)


def convert_strelka2_somatic_snv_vcf_to_dataframe(
        vcf_file: str,
        sequencing_platform: str,
        sample_id: str,
        tumor_sample_id: str,
        normal_sample_id: str
    ) -> pd.DataFrame:
    """
    Converts a Strelka2 (somatic-mode) SNV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   VCF file.
    sequencing_platform     :   Sequencing platform.
    sample_id               :   Sample ID.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID.

    Returns
    -------
    df                      :   DataFrame with the keys of
                                default_parameters.SMALL_VARIANT_ATTRIBUTES
                                as the columns
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['sample_id'] = sample_id
        curr_row['tumor_sample_id'] = tumor_sample_id
        curr_row['normal_sample_id'] = normal_sample_id
        curr_row['variant_id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.STRELKA2_SOMATIC
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['filter'] = curr_row['filter'].replace(';', ',')
        curr_row['quality_score'] = safely_retrieve_value(dict=row, key='QUAL', default_value=curr_row['quality_score'], type='float')


def convert_strelka2_germline_vcf_to_dataframe(
        vcf_file: str,
        sequencing_platform: str,
        sample_id: str,
        tumor_sample_id: str,
        normal_sample_id: str = ''
    ) -> pd.DataFrame:
    """
    Converts a Strelka2 VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   VCF file.
    sequencing_platform     :   Sequencing platform.
    sample_id               :   Sample ID.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID. If this parameter is
                                an empty string, it is assumed that the
                                variant calling was performed using a tumor
                                sample only.

    Returns
    -------
    df                      :   DataFrame with the keys of
                                default_parameters.SMALL_VARIANT_ATTRIBUTES
                                as the columns
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['sample_id'] = sample_id
        curr_row['tumor_sample_id'] = tumor_sample_id
        curr_row['normal_sample_id'] = normal_sample_id
        curr_row['variant_id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.STRELKA2_GERMLINE
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['filter'] = curr_row['filter'].replace(';', ',')
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

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'SNVHPOL':
                curr_row['snv_contextual_homopolymer_length'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'CIGAR':
                curr_row['cigar'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'RU':
                curr_row['smallest_repeating_sequence_unit'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'REFREP':
                curr_row['smallest_repeating_sequence_unit_reference_repeat_count'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'IDREP':
                curr_row['smallest_repeating_sequence_unit_allele_repeat_count'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'MQ':
                curr_row['mapping_quality_root_mean_square'] = float(curr_info_elements[1])

        # Extract FORMAT
        format = str(row['FORMAT']).split(':')
        tumor_sample = str(row[tumor_sample_id]).split(':')
        if normal_sample_id == '':
            normal_sample = ''
        else:
            normal_sample = str(row[normal_sample_id]).split(':')

        if 'GT' in format:
            curr_row['genotype'] = safely_retrieve_value(dict=tumor_sample, key=format.index('GT'), default_value=curr_row['genotype'], type='str')
            if normal_sample_id != '':
                curr_row['normal_genotype'] = safely_retrieve_value(dict=normal_sample, key=format.index('GT'), default_value=curr_row['normal_genotype'], type='str')
        if 'GQ' in format:
            curr_row['genotype_quality'] = safely_retrieve_value(dict=tumor_sample, key=format.index('GQ'), default_value=curr_row['genotype_quality'], type='int')
            if normal_sample_id != '':
                curr_row['normal_genotype_quality'] = safely_retrieve_value(dict=normal_sample, key=format.index('GQ'), default_value=curr_row['normal_genotype_quality'], type='int')
        if 'GQX' in format:
            curr_row['genotype_quality_recalibrated'] = safely_retrieve_value(dict=tumor_sample, key=format.index('GQX'), default_value=curr_row['genotype_quality_recalibrated'], type='int')
            if normal_sample_id != '':
                curr_row['normal_genotype_quality_recalibrated'] = safely_retrieve_value(dict=normal_sample, key=format.index('GQX'), default_value=curr_row['normal_genotype_quality_recalibrated'], type='int')
        if 'DP' in format:
            curr_row['total_depth'] = safely_retrieve_value(dict=tumor_sample, key=format.index('DP'), default_value=curr_row['total_depth'], type='int')
            if normal_sample_id != '':
                curr_row['normal_total_depth'] = safely_retrieve_value(dict=normal_sample, key=format.index('DP'), default_value=curr_row['normal_total_depth'], type='int')
        if 'DPF' in format:
            curr_row['filtered_basecalls_prior_to_genotyping'] = safely_retrieve_value(dict=tumor_sample, key=format.index('DPF'), default_value=curr_row['filtered_basecalls_prior_to_genotyping'], type='int')
            if normal_sample_id != '':
                curr_row['normal_filtered_basecalls_prior_to_genotyping'] = safely_retrieve_value(dict=normal_sample, key=format.index('DPF'), default_value=curr_row['normal_filtered_basecalls_prior_to_genotyping'], type='int')
        if 'AD' in format:
            curr_ad = safely_retrieve_value(dict=tumor_sample, key=format.index('AD'), default_value='unknown', type='str').split(',')
            try:
                curr_row['reference_reads_count'] = int(curr_ad[0])
            except:
                pass
            try:
                curr_row['variant_reads_count'] = int(curr_ad[1])
            except:
                pass
            if normal_sample_id != '':
                curr_ad_normal = safely_retrieve_value(dict=normal_sample, key=format.index('AD'), default_value='unknown', type='str').split(',')
                try:
                    curr_row['normal_reference_reads_count'] = int(curr_ad_normal[0])
                except:
                    pass
        if 'ADF' in format:
            curr_row['allelic_depths_forward_strand'] = safely_retrieve_value(dict=tumor_sample, key=format.index('ADF'), default_value=curr_row['allelic_depths_forward_strand'], type='int')
            if normal_sample_id != '':
                curr_row['normal_allelic_depths_forward_strand'] = safely_retrieve_value(dict=normal_sample, key=format.index('ADF'), default_value=curr_row['normal_allelic_depths_forward_strand'], type='int')
        if 'ADR' in format:
            curr_row['allelic_depths_reverse_strand'] = safely_retrieve_value(dict=tumor_sample, key=format.index('ADR'), default_value=curr_row['allelic_depths_reverse_strand'], type='int')
            if normal_sample_id != '':
                curr_row['normal_allelic_depths_reverse_strand'] = safely_retrieve_value(dict=normal_sample, key=format.index('ADR'), default_value=curr_row['normal_allelic_depths_reverse_strand'], type='int')
        if 'SB' in format:
            curr_row['strand_bias'] = safely_retrieve_value(dict=tumor_sample, key=format.index('SB'), default_value=curr_row['strand_bias'], type='float')
            if normal_sample_id != '':
                curr_row['normal_strand_bias'] = safely_retrieve_value(dict=normal_sample, key=format.index('SB'), default_value=curr_row['strand_bias'], type='float')
        if 'PL' in format:
            curr_row['phred_scale_genotype_likelihoods'] = safely_retrieve_value(dict=tumor_sample, key=format.index('PL'), default_value=curr_row['phred_scale_genotype_likelihoods'], type='str')
            if normal_sample_id != '':
                curr_row['normal_phred_scale_genotype_likelihoods'] = safely_retrieve_value(dict=normal_sample, key=format.index('PL'), default_value=curr_row['normal_phred_scale_genotype_likelihoods'], type='str')

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
            curr_row['variant_id'] = VariantCallingMethods.SmallVariantCallingMethods.STRELKA2_GERMLINE + '.' + \
                             curr_row['variant_type'] + '.' + \
                             str(curr_idx)
            curr_idx += 1

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df
