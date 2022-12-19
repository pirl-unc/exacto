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
handling VCF files.
"""


import gzip
import pandas as pd
from collections import defaultdict
from ..logging import get_logger
from ..constants import *
from ..default_parameters import *


logger = get_logger(__name__)


def __safely_retrieve_value(dict, key, default_value, type):
    """
    Safely converts a value from a VCF row.

    Parameters
    ----------
    dict            :   Dictionary.
    key             :   Key.
    default_value   :   Default value.
    type            :   Type ('str', 'int', 'float).

    Returns
    -------
    value   :   Value converted to 'type'.
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


def read_vcf_file(vcf_file: str) -> pd.DataFrame:
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


def convert_deepvariant_vcf_to_dataframe(vcf_file: str,
                                         sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a DeepVariant VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file            :   VCF file.
    sequencing_platform :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chrom'
    'pos'
    'ref'
    'alt'
    'filter'
    'quality_score'
    'variant_type'
    'variant_sequence'
    'variant_size'
    'genotype'
    'genotype_quality'
    'total_coverage'
    'reference_reads_count'
    'variant_reads_count'
    'variant_allele_fraction'
    'phred_scale_genotype_likelihoods'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = __safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = __safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = __safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = __safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = __safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['quality_score'] = __safely_retrieve_value(dict=row, key='QUAL', default_value=curr_row['quality_score'], type='float')

        # Make sure 'chr' is in chrom
        if 'chr' not in curr_row['chrom']:
            curr_row['chrom'] = 'chr' + curr_row['chrom']

        if len(curr_row['ref']) == 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SINGLE_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
            curr_row['variant_size'] = 1
        elif len(curr_row['ref']) == 1 and len(curr_row['alt']) > 1:
            if ',' in curr_row['alt']:
                curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
                curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
            else:
                curr_row['variant_type'] = SmallVariantTypes.SMALL_INSERTION
                curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
                curr_row['variant_size'] = len(curr_row['alt'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SMALL_DELETION
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='ref', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
            curr_row['variant_size'] = len(curr_row['ref'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) > 1:
            curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
        else:
            logger.warning('Unknown variant type. REF: %s. ALT: %s' %
                           (curr_row['ref'], curr_row['alt']))

        # Extract FORMAT
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['tumor_genotype'] = __safely_retrieve_value(dict=sample, key=format.index('GT'), default_value=curr_row['tumor_genotype'], type='str')
        if 'GQ' in format:
            curr_row['tumor_genotype_quality'] = __safely_retrieve_value(dict=sample, key=format.index('GQ'), default_value=curr_row['tumor_genotype_quality'], type='str')
        if 'DP' in format:
            curr_row['tumor_total_coverage'] = __safely_retrieve_value(dict=sample, key=format.index('DP'), default_value=curr_row['tumor_total_coverage'], type='int')
        if 'AD' in format:
            curr_ad = str(sample[format.index('AD')]).split(',')
            curr_row['tumor_reference_reads_count'] = int(curr_ad[0])
            if curr_row['variant_type'] == SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT:
                curr_row['tumor_variant_reads_count'] = ','.join(curr_ad[1:])
            else:
                curr_row['tumor_variant_reads_count'] = int(curr_ad[1])
        if 'VAF' in format:
            if curr_row['variant_type'] == SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT:
                curr_row['variant_allele_fraction'] = __safely_retrieve_value(dict=sample, key=format.index('VAF'), default_value=curr_row['variant_allele_fraction'], type='str')
            else:
                curr_row['variant_allele_fraction'] = __safely_retrieve_value(dict=sample, key=format.index('VAF'), default_value=curr_row['variant_allele_fraction'], type='float')
        if 'PL' in format:
            curr_row['phred_scale_genotype_likelihoods'] = __safely_retrieve_value(dict=sample, key=format.index('PL'), default_value=curr_row['phred_scale_genotype_likelihoods'], type='str')

        # Update total_coverage if it is currently unknown but can be inferred
        if type(curr_row['tumor_variant_reads_count']) == int and \
            type(curr_row['tumor_reference_reads_count']) == int and \
            curr_row['tumor_total_coverage'] == SMALL_VARIANT_ATTRIBUTES['tumor_total_coverage']:
            curr_row['tumor_total_coverage'] = curr_row['tumor_reference_reads_count'] + curr_row['tumor_variant_reads_count']

        # Update ID
        if curr_row['id'] == '.':
            curr_row['id'] = VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT + '.' + \
                             curr_row['variant_type'] + '.' + \
                             str(curr_idx)
            curr_idx += 1

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_gatk4_mutect2_vcf_to_dataframe(vcf_file: str,
                                           sequencing_platform: str,
                                           tumor_sample_id: str,
                                           normal_sample_id='') -> pd.DataFrame:
    """
    Convert a GATK4-Mutect2 VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   VCF file.
    sequencing_platform     :   Sequencing platform.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID. If this parameter is
                                an empty string, it is assumed that the
                                variant calling was performed using a tumor
                                sample only.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chrom'
    'pos'
    'ref'
    'alt'
    'variant_type'
    'filter'
    'as_sb_table'
    'ecnt'
    'germq'
    'mbq'
    'mfrl'
    'mmq'
    'mpos'
    'nalod'
    'nlod'
    'popaf'
    'tlod'
    'tumor_reference_reads_count'
    'tumor_variant_reads_count'
    'tumor_total_coverage'
    'variant_allele_fraction'
    'normal_reference_reads_count'
    'normal_total_coverage'
    'tumor_genotype'
    'normal_genotype'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = __safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = __safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = __safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = __safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = __safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['quality_score'] = __safely_retrieve_value(dict=row, key='QUAL', default_value=curr_row['quality_score'], type='float')

        # Make sure 'chr' is in chrom
        if 'chr' not in curr_row['chrom']:
            curr_row['chrom'] = 'chr' + curr_row['chrom']

        if len(curr_row['ref']) == 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SINGLE_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
            curr_row['variant_size'] = 1
        elif len(curr_row['ref']) == 1 and len(curr_row['alt']) > 1:
            if ',' in curr_row['alt']:
                curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
                curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
            else:
                curr_row['variant_type'] = SmallVariantTypes.SMALL_INSERTION
                curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
                curr_row['variant_size'] = len(curr_row['alt'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) == 1:
            curr_row['variant_type'] = SmallVariantTypes.SMALL_DELETION
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='ref', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
            curr_row['variant_size'] = len(curr_row['ref'][1:])
        elif len(curr_row['ref']) > 1 and len(curr_row['alt']) > 1:
            curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
            curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
        else:
            logger.warning('Unknown variant type. REF: %s. ALT: %s' %
                           (curr_row['ref'], curr_row['alt']))

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'AS_SB_TABLE':
                curr_row['allele_specific_strand_bias_table'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'ECNT':
                curr_row['haplotype_events'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'GERMQ':
                curr_row['alt_allele_germline_quality'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'MBQ':
                curr_row['allele_median_base_qualities'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'MFRL':
                curr_row['allele_median_fragment_length'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'MMQ':
                curr_row['allele_median_mapping_quality'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'MPOS':
                curr_row['median_distance_from_read_end'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'NALOD':
                curr_row['negative_log10_odds_artifact'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'NLOD':
                curr_row['log10_odds_artifact'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'POPAF':
                curr_row['negative_log_10_population'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'TLOD':
                curr_row['log10_likelihood_ratio_score_variant_exists'] = str(curr_info_elements[1])

        # Extract FORMAT
        format = str(row['FORMAT']).split(';')
        tumor_sample = str(row[tumor_sample_id]).split(';')
        if normal_sample_id == '':
            normal_sample = ''
        else:
            normal_sample = str(row[normal_sample_id]).split(';')

        if 'AD' in format:
            print(tumor_sample['AD'])
            curr_row['tumor_reference_reads_count'] = int(tumor_sample['AD'].split(',')[0])
            curr_row['tumor_variant_reads_count'] = int(tumor_sample['AD'].split(',')[1])
            if normal_sample_id != '':
                curr_row['normal_reference_reads_count'] = int(normal_sample['AD'].split(',')[0])
        if 'DP' in format:
            curr_row['tumor_total_coverage'] = int(tumor_sample['DP'])
            if normal_sample_id != '':
                curr_row['normal_total_coverage'] = int(normal_sample['DP'])
        if 'SB' in format:
            curr_row['tumor_strand_bias_fisher_exact_test_component_statistics'] = str(tumor_sample['SB'])
            if normal_sample_id != '':
                curr_row['normal_strand_bias_fisher_exact_test_component_statistics'] = str(normal_sample['SB'])
        if 'F1R2' in format:
            curr_row['tumor_f1r2_reads_count'] = str(tumor_sample['F1R2'])
            if normal_sample_id != '':
                curr_row['normal_f1r2_reads_count'] = str(normal_sample['F1R2'])
        if 'F2R1' in format:
            curr_row['tumor_f2r1_reads_count'] = str(tumor_sample['F2R1'])
            if normal_sample_id != '':
                curr_row['normal_f2r1_reads_count'] = str(normal_sample['F2R1'])
        if 'GT' in format:
            curr_row['tumor_genotype'] = str(tumor_sample['GT'])
            if normal_sample_id != '':
                curr_row['normal_genotype'] = str(normal_sample['GT'])

        # Update total_coverage if it is currently unknown but can be inferred
        if type(curr_row['tumor_variant_reads_count']) == int and \
            type(curr_row['tumor_reference_reads_count']) == int and \
            curr_row['tumor_total_coverage'] == SMALL_VARIANT_ATTRIBUTES['tumor_total_coverage']:
            curr_row['tumor_total_coverage'] = curr_row['tumor_reference_reads_count'] + curr_row['tumor_variant_reads_count']

        if type(curr_row['tumor_variant_reads_count']) == int and \
            type(curr_row['tumor_total_coverage']) == int:
            curr_row['variant_allele_fraction'] = float(curr_row['tumor_variant_reads_count']) / float(curr_row['tumor_total_coverage'])

        # Update ID
        if curr_row['id'] == '.':
            curr_row['id'] = VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2 + '.' + \
                             curr_row['variant_type'] + '.' + \
                             str(curr_idx)
            curr_idx += 1

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_strelka2_vcf_to_dataframe(vcf_file: str,
                                      sequencing_platform: str,
                                      tumor_sample_id: str,
                                      normal_sample_id='') -> pd.DataFrame:
    """
    Converts a Strelka2 VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   VCF file.
    sequencing_platform     :   Sequencing platform.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID. If this parameter is
                                an empty string, it is assumed that the
                                variant calling was performed using a tumor
                                sample only.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    list_data = []
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        curr_row = SMALL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.SmallVariantCallingMethods.STRELKA2
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chrom'] = __safely_retrieve_value(dict=row, key='CHROM', default_value=curr_row['chrom'], type='str')
        curr_row['pos'] = __safely_retrieve_value(dict=row, key='POS', default_value=curr_row['pos'], type='int')
        curr_row['ref'] = __safely_retrieve_value(dict=row, key='REF', default_value=curr_row['ref'], type='str').upper()
        curr_row['alt'] = __safely_retrieve_value(dict=row, key='ALT', default_value=curr_row['alt'], type='str').upper()
        curr_row['filter'] = __safely_retrieve_value(dict=row, key='FILTER', default_value=curr_row['filter'], type='str')
        curr_row['quality_score'] = __safely_retrieve_value(dict=row, key='QUAL', default_value=curr_row['quality_score'], type='float')

        # Make sure 'chr' is in chrom
        if 'chr' not in curr_row['chrom']:
            curr_row['chrom'] = 'chr' + curr_row['chrom']

    return pd.DataFrame()
    #     if len(curr_row['ref']) == 1 and len(curr_row['alt']) == 1:
    #         curr_row['variant_type'] = SmallVariantTypes.SINGLE_NUCLEOTIDE_VARIANT
    #         curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
    #         curr_row['variant_size'] = 1
    #     elif len(curr_row['ref']) == 1 and len(curr_row['alt']) > 1:
    #         if ',' in curr_row['alt']:
    #             curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
    #             curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
    #         else:
    #             curr_row['variant_type'] = SmallVariantTypes.SMALL_INSERTION
    #             curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
    #             curr_row['variant_size'] = len(curr_row['alt'][1:])
    #     elif len(curr_row['ref']) > 1 and len(curr_row['alt']) == 1:
    #         curr_row['variant_type'] = SmallVariantTypes.SMALL_DELETION
    #         curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='ref', default_value=curr_row['variant_sequence'], type='str')[1:].upper()
    #         curr_row['variant_size'] = len(curr_row['ref'][1:])
    #     elif len(curr_row['ref']) > 1 and len(curr_row['alt']) > 1:
    #         curr_row['variant_type'] = SmallVariantTypes.MULTI_NUCLEOTIDE_VARIANT
    #         curr_row['variant_sequence'] = __safely_retrieve_value(dict=row, key='alt', default_value=curr_row['variant_sequence'], type='str').upper()
    #     else:
    #         logger.warning('Unknown variant type. REF: %s. ALT: %s' %
    #                        (curr_row['ref'], curr_row['alt']))
    #
    #     # Extract INFO
    #     info = str(row['INFO']).split(';')
    #     for curr_info in info:
    #         curr_info_elements = curr_info.split('=')
    #         if curr_info_elements[0] == 'AS_SB_TABLE':
    #             curr_row['allele_specific_strand_bias_table'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'ECNT':
    #             curr_row['haplotype_events'] = int(curr_info_elements[1])
    #         if curr_info_elements[0] == 'GERMQ':
    #             curr_row['alt_allele_germline_quality'] = int(curr_info_elements[1])
    #         if curr_info_elements[0] == 'MBQ':
    #             curr_row['allele_median_base_qualities'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'MFRL':
    #             curr_row['allele_median_fragment_length'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'MMQ':
    #             curr_row['allele_median_mapping_quality'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'MPOS':
    #             curr_row['median_distance_from_read_end'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'NALOD':
    #             curr_row['negative_log10_odds_artifact'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'NLOD':
    #             curr_row['log10_odds_artifact'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'POPAF':
    #             curr_row['negative_log_10_population'] = str(curr_info_elements[1])
    #         if curr_info_elements[0] == 'TLOD':
    #             curr_row['log10_likelihood_ratio_score_variant_exists'] = str(curr_info_elements[1])
    #
    #     # Extract FORMAT
    #     format = str(row['FORMAT']).split(';')
    #     tumor_sample = str(row[tumor_sample_id]).split(';')
    #     if normal_sample_id == '':
    #         normal_sample = ''
    #     else:
    #         normal_sample = str(row[normal_sample_id]).split(';')
    #
    #     if 'AD' in format:
    #         print(tumor_sample['AD'])
    #         curr_row['tumor_reference_reads_count'] = int(tumor_sample['AD'].split(',')[0])
    #         curr_row['tumor_variant_reads_count'] = int(tumor_sample['AD'].split(',')[1])
    #         if normal_sample_id != '':
    #             curr_row['normal_reference_reads_count'] = int(normal_sample['AD'].split(',')[0])
    #     if 'DP' in format:
    #         curr_row['tumor_total_coverage'] = int(tumor_sample['DP'])
    #         if normal_sample_id != '':
    #             curr_row['normal_total_coverage'] = int(normal_sample['DP'])
    #     if 'SB' in format:
    #         curr_row['tumor_strand_bias_fisher_exact_test_component_statistics'] = str(tumor_sample['SB'])
    #         if normal_sample_id != '':
    #             curr_row['normal_strand_bias_fisher_exact_test_component_statistics'] = str(normal_sample['SB'])
    #     if 'F1R2' in format:
    #         curr_row['tumor_f1r2_reads_count'] = str(tumor_sample['F1R2'])
    #         if normal_sample_id != '':
    #             curr_row['normal_f1r2_reads_count'] = str(normal_sample['F1R2'])
    #     if 'F2R1' in format:
    #         curr_row['tumor_f2r1_reads_count'] = str(tumor_sample['F2R1'])
    #         if normal_sample_id != '':
    #             curr_row['normal_f2r1_reads_count'] = str(normal_sample['F2R1'])
    #     if 'GT' in format:
    #         curr_row['tumor_genotype'] = str(tumor_sample['GT'])
    #         if normal_sample_id != '':
    #             curr_row['normal_genotype'] = str(normal_sample['GT'])
    #
    #     # Update total_coverage if it is currently unknown but can be inferred
    #     if type(curr_row['tumor_variant_reads_count']) == int and \
    #         type(curr_row['tumor_reference_reads_count']) == int and \
    #         curr_row['tumor_total_coverage'] == SMALL_VARIANT_ATTRIBUTES['tumor_total_coverage']:
    #         curr_row['tumor_total_coverage'] = curr_row['tumor_reference_reads_count'] + curr_row['tumor_variant_reads_count']
    #
    #     if type(curr_row['tumor_variant_reads_count']) == int and \
    #         type(curr_row['tumor_total_coverage']) == int:
    #         curr_row['variant_allele_fraction'] = float(curr_row['tumor_variant_reads_count']) / float(curr_row['tumor_total_coverage'])
    #
    #     # Update ID
    #     if curr_row['id'] == '.':
    #         curr_row['id'] = VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2 + '.' + \
    #                          curr_row['variant_type'] + '.' + \
    #                          str(curr_idx)
    #         curr_idx += 1
    #
    #     # Append to list
    #     list_data.append(curr_row)
    #
    # df = pd.DataFrame.from_dict(list_data)
    # logger.info('%i rows in the returning DataFrame.' % len(df))
    # return df


def convert_sniffles2_vcf_to_dataframe(vcf_file: str,
                                       sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a Sniffles2 VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    for row in df_vcf.to_dict('records'):
        curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chr_1'] = str(row['CHROM'])
        curr_row['chr_2'] = str(row['CHROM'])
        curr_row['pos_1'] = int(row['POS'])
        curr_row['ref'] = str(row['REF']).upper()
        curr_row['alt'] = str(row['ALT']).upper()
        curr_row['filter'] = str(row['FILTER'])
        if row['QUAL'] != '.':
            curr_row['quality_score'] = float(row['QUAL'])

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'PRECISE':
                curr_row['is_precise'] = True
            if curr_info_elements[0] == 'IMPRECISE':
                curr_row['is_precise'] = False
            if curr_info_elements[0] == 'SVTYPE':
                curr_row['sv_type'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'SVLEN':
                curr_row['sv_size'] = abs(int(curr_info_elements[1]))
            if curr_info_elements[0] == 'END':
                curr_row['pos_2'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'RNAMES':
                curr_row['read_ids'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'COVERAGE':
                curr_row['coverage'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'STRAND':
                curr_row['strand'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'AF':
                curr_row['variant_allele_fraction'] = float(curr_info_elements[1])
            if curr_info_elements[0] == 'STDEV_LEN':
                curr_row['sv_size_stdev'] = float(curr_info_elements[1])
            if curr_info_elements[0] == 'STDEV_POS':
                curr_row['sv_pos_stdev'] = float(curr_info_elements[1])
            if curr_info_elements[0] == 'CHR2':
                curr_row['chr_2'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'NM':
                curr_row['nm'] = float(curr_info_elements[1])
            if curr_info_elements[0] == 'SUPPORT_LONG':
                curr_row['support_long'] = int(curr_info_elements[1])

        # Make sure 'chr' is in chr_1 and chr_2
        if 'chr' not in curr_row['chr_1'] and curr_row['chr_1'] != '':
            curr_row['chr_1'] = 'chr' + curr_row['chr_1']
        if 'chr' not in curr_row['chr_2'] and curr_row['chr_2'] != '':
            curr_row['chr_2'] = 'chr' + curr_row['chr_2']

        # Update position 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            alt_val = curr_row['alt'].split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            curr_row['pos_2'] = int(alt_val)

        # Update SV size for 'BND'
        if curr_row['sv_type'] == 'BND' and curr_row['chr_1'] == curr_row['chr_2']:
            curr_row['sv_size'] = abs(curr_row['pos_1'] - curr_row['pos_2'])

        # Update insertion sequence
        if curr_row['sv_type'] == 'INS':
            curr_row['insertion_sequence'] = curr_row['alt']

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = str(sample[format.index('GT')])
        if 'GQ' in format:
            curr_row['genotype_quality'] = float(sample[format.index('GQ')])
        if 'DR' in format:
            curr_row['reference_reads_count'] = int(sample[format.index('DR')])
        if 'DV' in format:
            curr_row['variant_reads_count'] = int(sample[format.index('DV')])
        if curr_row['variant_reads_count'] > 0 and curr_row['reference_reads_count'] >= 0:
            curr_row['total_coverage'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_cutesv_vcf_to_dataframe(vcf_file: str,
                                    sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a CuteSV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    for row in df_vcf.to_dict('records'):
        curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.StructuralVariantCallingMethods.CUTESV
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chr_1'] = str(row['CHROM'])
        curr_row['chr_2'] = str(row['CHROM'])
        curr_row['pos_1'] = int(row['POS'])
        curr_row['ref'] = str(row['REF']).upper()
        curr_row['alt'] = str(row['ALT']).upper()
        curr_row['filter'] = str(row['FILTER'])
        if row['QUAL'] != '.':
            curr_row['quality_score'] = float(row['QUAL'])

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'PRECISE':
                curr_row['is_precise'] = True
            if curr_info_elements[0] == 'IMPRECISE':
                curr_row['is_precise'] = False
            if curr_info_elements[0] == 'SVTYPE':
                try:
                    curr_row['sv_type'] = str(curr_info_elements[1])
                except:
                    pass
            if curr_info_elements[0] == 'SVLEN':
                curr_row['sv_size'] = abs(int(curr_info_elements[1]))
            if curr_info_elements[0] == 'RE':
                curr_row['variant_reads_count'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'RNAMES':
                curr_row['read_ids'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'END':
                curr_row['pos_2'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'STRAND':
                curr_row['strand'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'CIPOS':
                curr_row['ci_pos'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'CILEN':
                curr_row['ci_len'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'AF' and curr_info_elements[0] != '.':
                curr_row['variant_allele_fraction'] = float(curr_info_elements[1])

        # Update chromosome 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            alt_val = curr_row['alt'].split(":")[0]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            curr_row['chr_2'] = str(alt_val)

        # Update position 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            alt_val = curr_row['alt'].split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            curr_row['pos_2'] = int(alt_val)

        # Update SV size for 'BND'
        if curr_row['sv_type'] == 'BND' and curr_row['chr_1'] == curr_row['chr_2']:
            curr_row['sv_size'] = abs(curr_row['pos_2'] - curr_row['pos_1'])

        # Make sure 'chr' is in chr_1 and chr_2
        if 'chr' not in curr_row['chr_1'] and curr_row['chr_1'] != '':
            curr_row['chr_1'] = 'chr' + curr_row['chr_1']
        if 'chr' not in curr_row['chr_2'] and curr_row['chr_2'] != '':
            curr_row['chr_2'] = 'chr' + curr_row['chr_2']

        # Update strand for 'BND'
        if curr_row['sv_type'] == 'BND':
            if curr_row['alt'][0:2] == 'N[':
                curr_row['strand'] = '+-'
            elif curr_row['alt'][0:2] == 'N]':
                curr_row['strand'] = '-+'
            elif curr_row['alt'][-2:] == ']N':
                curr_row['strand'] = '+-'
            elif curr_row['alt'][-2:] == '[N':
                curr_row['strand'] = '-+'
            else:
                curr_row['strand'] = ''

        # Update insertion sequence
        if curr_row['sv_type'] == StructuralVariantTypes.INSERTION:
            curr_row['insertion_sequence'] = curr_row['alt'][1:]

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = str(sample[format.index('GT')])
        if 'GQ' in format:
            curr_row['genotype_quality'] = float(sample[format.index('GQ')])
        if 'DR' in format:
            curr_row['reference_reads_count'] = int(sample[format.index('DR')])
        if 'DV' in format:
            curr_row['variant_reads_count'] = int(sample[format.index('DV')])
        if curr_row['variant_reads_count'] > 0 and curr_row['reference_reads_count'] >= 0:
            curr_row['total_coverage'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_svim_vcf_to_dataframe(vcf_file: str,
                                  sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a SVIM VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'sv_type'
    'sv_size'
    'sv_pos_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    for row in df_vcf.to_dict('records'):
        curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
        curr_row['variant_calling_method'] = VariantCallingMethods.StructuralVariantCallingMethods.SVIM
        curr_row['sequencing_platform'] = sequencing_platform
        curr_row['chr_1'] = str(row['CHROM'])
        curr_row['chr_2'] = str(row['CHROM'])
        curr_row['pos_1'] = int(row['POS'])
        curr_row['ref'] = str(row['REF']).upper()
        curr_row['alt'] = str(row['ALT']).upper()
        curr_row['filter'] = str(row['FILTER'])
        curr_row['is_precise'] = 'unknown'
        if row['QUAL'] != '.':
            curr_row['quality_score'] = float(row['QUAL'])

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            curr_info_elements = curr_info.split('=')
            if curr_info_elements[0] == 'SVTYPE':
                curr_row['sv_type'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'SVLEN':
                curr_row['sv_size'] = abs(int(curr_info_elements[1]))
            if curr_info_elements[0] == 'END':
                curr_row['pos_2'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'SUPPORT':
                curr_row['variant_reads_count'] = int(curr_info_elements[1])
            if curr_info_elements[0] == 'READS':
                curr_row['read_ids'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'SEQS':
                curr_row['insertion_sequence'] = str(curr_info_elements[1])
            if curr_info_elements[0] == 'STD_SPAN':
                try:
                    curr_row['std_span'] = float(curr_info_elements[1])
                except:
                    curr_row['std_span'] = ''
            if curr_info_elements[0] == 'STD_POS':
                try:
                    curr_row['sv_pos_stdev'] = float(curr_info_elements[1])
                except:
                    curr_row['sv_pos_stdev'] = ''

        # Convert 'INVDUP' to 'DUP'
        if "DUP" in curr_row['sv_type']:
            curr_row['sv_type'] = 'DUP'

        # Update chromosome 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            alt_val = curr_row['alt'].split(":")[0]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            curr_row['chr_2'] = str(alt_val)

        # Update position 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            alt_val = curr_row['alt'].split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            curr_row['pos_2'] = int(alt_val)

        # Update SV size for 'BND'
        if curr_row['sv_type'] == 'BND' and curr_row['chr_1'] == curr_row['chr_2']:
            curr_row['sv_size'] = abs(curr_row['pos_2'] - curr_row['pos_1'])

        # Update strand for 'BND'
        if curr_row['sv_type'] == 'BND':
            if curr_row['alt'][0:2] == 'N[':
                curr_row['strand'] = '+-'
            elif curr_row['alt'][0:2] == 'N]':
                curr_row['strand'] = '-+'
            elif curr_row['alt'][-2:] == ']N':
                curr_row['strand'] = '+-'
            elif curr_row['alt'][-2:] == '[N':
                curr_row['strand'] = '-+'
            else:
                curr_row['strand'] = ''

        # Update insertion sequence
        if curr_row['sv_type'] == 'INS':
            curr_row['insertion_sequence'] = curr_row['alt'][1:]

        # Make sure 'chr' is in chr_1 and chr_2
        if 'chr' not in curr_row['chr_1'] and curr_row['chr_1'] != '':
            curr_row['chr_1'] = 'chr' + curr_row['chr_1']
        if 'chr' not in curr_row['chr_2'] and curr_row['chr_2'] != '':
            curr_row['chr_2'] = 'chr' + curr_row['chr_2']

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = str(sample[format.index('GT')])
        if 'DP' in format:
            try:
                curr_row['total_coverage'] = int(sample[format.index('DP')])
            except:
                curr_row['total_coverage'] = -1
                logger.warning('Total coverage is not present for %s. Saving total coverage as -1.' % curr_row['id'])
        if 'AD' in format:
            try:
                curr_row['variant_reads_count'] = int(sample[format.index('AD')].split(',')[1])
            except:
                curr_row['variant_reads_count'] = -1
                logger.warning('Variant reads count is not present for %s. Saving variant reads count as -1.' % curr_row['id'])
            try:
                curr_row['reference_reads_count'] = int(sample[format.index('AD')].split(',')[0])
            except:
                curr_row['reference_reads_count'] = -1
                logger.warning('Reference reads count is not present for %s. Saving reference reads count as -1.' % curr_row['id'])
        if 'CN' in format:
            curr_row['tandem_duplication_copy_number'] = int(sample[format.index('CN')])
        if curr_row['variant_reads_count'] > 0 and curr_row['reference_reads_count'] >= 0:
            curr_row['total_coverage'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']
            curr_row['variant_allele_fraction'] = curr_row['variant_reads_count'] / curr_row['total_coverage']

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_pbsv_vcf_to_dataframe(vcf_file: str,
                                  sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a PBSV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
    """
    df_vcf = read_vcf_file(vcf_file=vcf_file)
    sample_key = df_vcf.columns.values.tolist()[-1]
    list_data = []
    included_mate_ids = set()
    for row in df_vcf.to_dict('records'):
        curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
        curr_row['id'] = str(row['ID'])
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

        if curr_row['sv_type'] == 'BND':
            # Check if current ID has been included
            if mate_id in included_mate_ids:
                continue
            included_mate_ids.add(curr_row['id'])

        # Update chromosome 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            curr_id = curr_row['id'].split("-")[1]
            curr_id = curr_id.split(":")[0]
            curr_row['chr_2'] = str(curr_id)

        # Update position 2 for 'BND'
        if curr_row['sv_type'] == 'BND':
            curr_id = curr_row['id'].split("-")[1]
            curr_id = curr_id.split(":")[1]
            curr_row['pos_2'] = int(curr_id)

        # Make sure 'chr' is in chr_1 and chr_2
        if 'chr' not in curr_row['chr_1'] and curr_row['chr_1'] != '':
            curr_row['chr_1'] = 'chr' + curr_row['chr_1']
        if 'chr' not in curr_row['chr_2'] and curr_row['chr_2'] != '':
            curr_row['chr_2'] = 'chr' + curr_row['chr_2']

        # Update SV size for 'BND'
        if curr_row['sv_type'] == 'BND' and curr_row['chr_1'] == curr_row['chr_2']:
            curr_row['sv_size'] = abs(curr_row['pos_2'] - curr_row['pos_1'])

        # Update strand for 'BND'
        if curr_row['sv_type'] == 'BND':
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
        if curr_row['sv_type'] == 'INS':
            curr_row['insertion_sequence'] = curr_row['alt']

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        sample = str(row[sample_key]).split(':')
        if 'GT' in format:
            curr_row['genotype'] = str(sample[format.index('GT')])
        if 'DP' in format:
            try:
                curr_row['total_coverage'] = int(sample[format.index('DP')])
            except:
                curr_row['total_coverage'] = -1
                logger.warning('Total coverage is not present for %s. Saving total coverage as -1.' % curr_row['id'])
        if 'AD' in format:
            try:
                curr_row['variant_reads_count'] = int(sample[format.index('AD')].split(',')[1])
            except:
                curr_row['variant_reads_count'] = -1
                logger.warning(
                    'Variant reads count is not present for %s. Saving variant reads count as -1.' % curr_row['id'])
            try:
                curr_row['reference_reads_count'] = int(sample[format.index('AD')].split(',')[0])
            except:
                curr_row['reference_reads_count'] = -1
                logger.warning(
                    'Reference reads count is not present for %s. Saving reference reads count as -1.' % curr_row['id'])
        if 'SAC' in format:
            curr_row['strand_reads'] = str(sample[format.index('SAC')])

        if curr_row['variant_reads_count'] > 0 and curr_row['reference_reads_count'] >= 0:
            curr_row['total_coverage'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']
            curr_row['variant_allele_fraction'] = curr_row['variant_reads_count'] / curr_row['total_coverage']

        # Append to list
        list_data.append(curr_row)

    df = pd.DataFrame.from_dict(list_data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df

