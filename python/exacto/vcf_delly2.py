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


# """
# The purpose of this python3 script is to implement functions
# related to handling DELLY2 VCF files.
# """
#
#
# import gzip
# from collections import defaultdict
# from .logging import get_logger
# from .constants import *
# from .default_parameters import *
# from .vcf import *
# from .variant import Variant
# from .variant_call import VariantCall
# from .variants_list import VariantsList
#
#
# logger = get_logger(__name__)
#
#
# def parse_delly2_callset(
#         df_vcf: pd.DataFrame,
#         sequencing_platform: str,
#         source_id: str
#     ) -> pd.DataFrame:
#     """
#     Parses a Delly2 DataFrame and
#     returns an instance of the VariantsList class.
#
#     Parameters
#     ----------
#     df_vcf                  :   DataFrame of rows from a Delly2 VCF file.
#     sequencing_platform     :   Sequencing platform.
#     source_id               :   Source ID.
#
#     Returns
#     -------
#     variants_list       :   An instance of the VariantsList class.
#     """
#     variants_list = VariantsList()
#     tumor_sample_id = df_vcf.columns.values.tolist()[-1]
#     curr_idx = 1
#     for row in df_vcf.to_dict('records'):
#         chr_1 = safely_retrieve_value(dict=row, key='CHROM', default_value=DEFAULT_ATTRIBUTE_VALUE, type=str)
#         pos_1 = safely_retrieve_value(dict=row, key='POS', default_value=DEFAULT_ATTRIBUTE_VALUE, type=int)
#         chr_2 = safely_retrieve_value(dict=row, key='CHROM', default_value=DEFAULT_ATTRIBUTE_VALUE, type=str)
#         ref = safely_retrieve_value(dict=row, key='REF', default_value=DEFAULT_ATTRIBUTE_VALUE, type=str)
#         alt = safely_retrieve_value(dict=row, key='ALT', default_value=DEFAULT_ATTRIBUTE_VALUE, type=str)
#         filter = safely_retrieve_value(dict=row, key='FILTER', default_value=DEFAULT_ATTRIBUTE_VALUE, type=str)
#         quality_score = safely_retrieve_value(dict=row, key='QUAL', default_value=DEFAULT_ATTRIBUTE_VALUE, type=float)
#         precise = False
#
# def convert_delly2_vcf_to_dataframe(
#         vcf_file: str,
#         sequencing_platform: str,
#         sample_id: str,
#         tumor_sample_id: str,
#         normal_sample_id: str = ''
#     ) -> pd.DataFrame:
#     """
#     Convert a Sniffles2 VCF file to a DataFrame.
#
#     Parameters
#     ----------
#     vcf_file                :   VCF file.
#     sequencing_platform     :   Sequencing platform.
#     sample_id               :   Sample ID.
#     tumor_sample_id         :   Tumor sample ID.
#     normal_sample_id        :   Normal sample ID. If this parameter is
#                                 an empty string, it is assumed that the
#                                 variant calling was performed using a tumor
#                                 sample only.
#
#     Returns
#     -------
#     df                      :   DataFrame with the keys of
#                                 default_parameters.STRUCTURAL_VARIANT_ATTRIBUTES
#                                 as the columns
#     """
#     df_vcf = read_vcf_file(vcf_file=vcf_file)
#     sample_key = df_vcf.columns.values.tolist()[-1]
#     list_data = []
#     curr_idx = 1
#     for row in df_vcf.to_dict('records'):
#         curr_row = STRUCTURAL_VARIANT_ATTRIBUTES.copy()
#         curr_row['sample_id'] = sample_id
#         curr_row['variant_id'] = str(row['ID'])
#         curr_row['variant_calling_method'] = VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2
#         curr_row['sequencing_platform'] = sequencing_platform
#         curr_row['chr_1'] = str(row['CHROM'])
#         curr_row['chr_2'] = str(row['CHROM'])
#         curr_row['pos_1'] = int(row['POS'])
#         curr_row['ref'] = str(row['REF']).upper()
#         curr_row['alt'] = str(row['ALT']).upper()
#         curr_row['filter'] = str(row['FILTER'])
#         if row['QUAL'] != '.':
#             curr_row['quality_score'] = float(row['QUAL'])
#
#         # Extract INFO
#         info = str(row['INFO']).split(';')
#         for curr_info in info:
#             curr_info_elements = curr_info.split('=')
#             if curr_info_elements[0] == 'PRECISE':
#                 curr_row['is_precise'] = True
#             if curr_info_elements[0] == 'IMPRECISE':
#                 curr_row['is_precise'] = False
#             if curr_info_elements[0] == 'SVTYPE':
#                 curr_row['sv_type'] = str(curr_info_elements[1])
#             if curr_info_elements[0] == 'SVLEN':
#                 curr_row['sv_size'] = abs(int(curr_info_elements[1]))
#             if curr_info_elements[0] == 'END':
#                 curr_row['pos_2'] = int(curr_info_elements[1])
#             if curr_info_elements[0] == 'RNAMES':
#                 curr_row['read_ids'] = str(curr_info_elements[1])
#             if curr_info_elements[0] == 'COVERAGE':
#                 curr_row['coverage'] = str(curr_info_elements[1])
#             if curr_info_elements[0] == 'STRAND':
#                 curr_row['strand'] = str(curr_info_elements[1])
#             if curr_info_elements[0] == 'AF':
#                 curr_row['variant_allele_fraction'] = float(curr_info_elements[1])
#             if curr_info_elements[0] == 'STDEV_LEN':
#                 curr_row['sv_size_stdev'] = float(curr_info_elements[1])
#             if curr_info_elements[0] == 'STDEV_POS':
#                 curr_row['sv_pos_stdev'] = float(curr_info_elements[1])
#             if curr_info_elements[0] == 'CHR2':
#                 curr_row['chr_2'] = str(curr_info_elements[1])
#             if curr_info_elements[0] == 'NM':
#                 curr_row['nm'] = float(curr_info_elements[1])
#             if curr_info_elements[0] == 'SUPPORT_LONG':
#                 curr_row['support_long'] = int(curr_info_elements[1])
#
#         # Make sure 'chr' is in chr_1 and chr_2
#         if ('chr' not in curr_row['chr_1']) and ('CHR' not in curr_row['chr_1']) and (curr_row['chr_1'] != ''):
#             curr_row['chr_1'] = 'chr' + curr_row['chr_1']
#         if ('chr' not in curr_row['chr_2']) and ('CHR' not in curr_row['chr_2']) and (curr_row['chr_2'] != ''):
#             curr_row['chr_2'] = 'chr' + curr_row['chr_2']
#
#         # Update position 2 for 'BND'
#         if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT:
#             alt_val = str(row['ALT']).split(":")[1]
#             alt_val = alt_val.replace("[", "")
#             alt_val = alt_val.replace("]", "")
#             alt_val = alt_val.replace("N", "")
#             curr_row['pos_2'] = int(alt_val)
#
#         # Update SV size for 'BND'
#         if curr_row['sv_type'] == StructuralVariantTypes.BREAKPOINT and curr_row['chr_1'] == curr_row['chr_2']:
#             curr_row['sv_size'] = abs(curr_row['pos_1'] - curr_row['pos_2'])
#
#         # Update insertion sequence
#         if curr_row['sv_type'] == StructuralVariantTypes.INSERTION:
#             curr_row['insertion_sequence'] = curr_row['alt']
#         else:
#             curr_row['insertion_sequence'] = ''
#
#         # Extract FORMAT (sample)
#         format = str(row['FORMAT']).split(':')
#         sample = str(row[sample_key]).split(':')
#         if 'GT' in format:
#             curr_row['genotype'] = safely_retrieve_value(dict=sample, key=format.index('GT'), default_value=curr_row['genotype'], type='str')
#         if 'GQ' in format:
#             curr_row['genotype_quality'] = safely_retrieve_value(dict=sample, key=format.index('GQ'), default_value=curr_row['genotype_quality'], type='float')
#         if 'DR' in format:
#             curr_row['reference_reads_count'] = safely_retrieve_value(dict=sample, key=format.index('DR'), default_value=curr_row['reference_reads_count'], type='int')
#         if 'DV' in format:
#             curr_row['variant_reads_count'] = safely_retrieve_value(dict=sample, key=format.index('DV'), default_value=curr_row['variant_reads_count'], type='int')
#
#         # Update total_depth if it is currently unknown but can be inferred
#         if type(curr_row['variant_reads_count']) == int and \
#             type(curr_row['reference_reads_count']) == int and \
#             curr_row['total_depth'] == STRUCTURAL_VARIANT_ATTRIBUTES['total_depth']:
#             curr_row['total_depth'] = curr_row['reference_reads_count'] + curr_row['variant_reads_count']
#
#         # Update variant_allele_fraction if it is currently unknown
#         if type(curr_row['variant_reads_count']) == int and \
#             type(curr_row['total_depth']) == int and \
#             curr_row['variant_allele_fraction'] == STRUCTURAL_VARIANT_ATTRIBUTES['variant_allele_fraction']:
#             if curr_row['total_depth'] > 0:
#                 curr_row['variant_allele_fraction'] = float(curr_row['variant_reads_count']) / float(curr_row['total_depth'])
#
#         # Update ID
#         if curr_row['variant_id'] == '.':
#             curr_row['variant_id'] = VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2 + '.' + \
#                              curr_row['variant_type'] + '.' + \
#                              str(curr_idx)
#             curr_idx += 1
#
#         # Append to list
#         list_data.append(curr_row)
#
#     df = pd.DataFrame.from_dict(list_data)
#     logger.info('%i rows in the returning DataFrame.' % len(df))
#     return df
#
