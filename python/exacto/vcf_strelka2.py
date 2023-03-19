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
from .logging import get_logger
from .constants import *
from .default_parameters import *
from .common import safely_retrieve_value, safely_convert_value
from .variant import Variant
from .variant_call import VariantCall
from .variants_list import VariantsList


logger = get_logger(__name__)


def parse_strelka2_somatic_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str,
        tumor_sample_id: str,
        normal_sample_id: str
    ) -> VariantsList:
    """
    Parses a Strelka2-Somatic DataFrame and
    returns an instance of the VariantsList class.

    Parameters
    ----------
    df_vcf                  :   DataFrame of rows from a Strelka2 (somatic mode) VCF file.
    sequencing_platform     :   Sequencing platform.
    source_id               :   Source ID.
    tumor_sample_id         :   Tumor sample ID.
    normal_sample_id        :   Normal sample ID.

    Returns
    -------
    variants_list           :   An instance of the VariantsList class.
    """
    variants_list = VariantsList()
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        variant_call = VariantCall()
        variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_id'] = safely_retrieve_value(
            dict=row,
            key='ID',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.nucleic_acid = NucleicAcidTypes.DNA
        variant_call.source_id = source_id
        variant_call.tumor_sample_id = tumor_sample_id
        variant_call.normal_sample_id = normal_sample_id
        variant_call.variant_calling_method = VariantCallingMethods.STRELKA2_SOMATIC
        variant_call.sequencing_platform = sequencing_platform
        variant_call.chr_1 = safely_retrieve_value(
            dict=row,
            key='CHROM',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.pos_1 = safely_retrieve_value(
            dict=row,
            key='POS',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=int
        )
        variant_call.chr_2 = safely_retrieve_value(
            dict=row,
            key='CHROM',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.ref = safely_retrieve_value(
            dict=row,
            key='REF',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.alt = safely_retrieve_value(
            dict=row,
            key='ALT',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.filter = safely_retrieve_value(
            dict=row,
            key='FILTER',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.quality_score = safely_retrieve_value(
            dict=row,
            key='QUAL',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=float
        )

        if len(variant_call.ref) == 1 and len(variant_call.alt) == 1:
            variant_call.variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
            variant_call.variant_sequences.append(variant_call.alt)
            variant_call.variant_size = 1
            variant_call.pos_2 = variant_call.pos_1
        elif len(variant_call.ref) == 1 and len(variant_call.alt) > 1:
            variant_call.variant_type = VariantTypes.INSERTION
            variant_call.variant_sequences.append(variant_call.alt[1:])
            variant_call.variant_size = len(variant_call.alt[1:])
            variant_call.pos_2 = variant_call.pos_1
        elif len(variant_call.ref) > 1 and len(variant_call.alt) == 1:
            variant_call.variant_type = VariantTypes.DELETION
            variant_call.variant_sequences.append(variant_call.ref[1:])
            variant_call.variant_size = len(variant_call.ref[1:])
            variant_call.pos_1 = variant_call.pos_1 + 1
            variant_call.pos_2 = variant_call.pos_1 + variant_call.variant_size - 1
        elif len(variant_call.ref) > 1 and len(variant_call.alt) > 1:
            variant_call.variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
            variant_call.variant_sequences.append(variant_call.alt)
            variant_call.variant_size = len(variant_call.alt)
            variant_call.pos_2 = variant_call.pos_1 + len(variant_call.variant_sequences[0]) - 1
        else:
            logger.error('Unknown variant type. REF: %s. ALT: %s' % (variant_call.ref, variant_call.alt))
            exit(1)

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            if '=' in curr_info:
                curr_info_elements = curr_info.split('=')
                variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_' + curr_info_elements[0].lower()] = safely_convert_value(
                    value=curr_info_elements[1],
                    default_value=DEFAULT_ATTRIBUTE_VALUE,
                    type=VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[VariantCallingMethods.STRELKA2_SOMATIC + '_' + curr_info_elements[0].lower()]
                )
            else:
                variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_' + curr_info.lower()] = True

        # Extract FORMAT
        format = str(row['FORMAT']).split(':')
        tumor_sample = str(row['TUMOR']).split(':')
        normal_sample = str(row['NORMAL']).split(':')
        for curr_format in format:
            variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_tumor_' + curr_format.lower()] = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(curr_format),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[VariantCallingMethods.STRELKA2_SOMATIC + '_' + curr_format.lower()]
            )
            variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_normal_' + curr_format.lower()] = safely_retrieve_value(
                dict=normal_sample,
                key=format.index(curr_format),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[VariantCallingMethods.STRELKA2_SOMATIC + '_' + curr_format.lower()]
            )

        # Update total_tumor_reads
        if VariantCallingMethods.STRELKA2_SOMATIC + '_tumor_dp' in variant_call.tool_attributes.keys():
            variant_call.total_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_tumor_dp']

        # Update total_normal_reads
        if VariantCallingMethods.STRELKA2_SOMATIC + '_normal_dp' in variant_call.tool_attributes.keys():
            variant_call.total_normal_reads = variant_call.tool_attributes[VariantCallingMethods.STRELKA2_SOMATIC + '_normal_dp']

        # Update ref_tumor_reads, alt_tumor_reads ref_normal_reads, alt_normal_reads
        if variant_call.variant_type == VariantTypes.SINGLE_NUCLEOTIDE_VARIANT:
            ref_format_value = ''
            if variant_call.ref == 'A':
                ref_format_value = 'AU'
            if variant_call.ref == 'C':
                ref_format_value = 'CU'
            if variant_call.ref == 'G':
                ref_format_value = 'GU'
            if variant_call.ref == 'T':
                ref_format_value = 'TU'

            alt_format_value = ''
            if variant_call.alt == 'A':
                alt_format_value = 'AU'
            if variant_call.alt == 'C':
                alt_format_value = 'CU'
            if variant_call.alt == 'G':
                alt_format_value = 'GU'
            if variant_call.alt == 'T':
                alt_format_value = 'TU'

            ref_tumor_reads = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(ref_format_value),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            ref_normal_reads = safely_retrieve_value(
                dict=normal_sample,
                key=format.index(ref_format_value),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            alt_tumor_reads = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(alt_format_value),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            alt_normal_reads = safely_retrieve_value(
                dict=normal_sample,
                key=format.index(alt_format_value),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            variant_call.ref_tumor_reads = safely_convert_value(
                value=ref_tumor_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
            variant_call.ref_normal_reads = safely_convert_value(
                value=ref_normal_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
            variant_call.alt_tumor_reads = safely_convert_value(
                value=alt_tumor_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
            variant_call.alt_normal_reads = safely_convert_value(
                value=alt_normal_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
        if variant_call.variant_type == VariantTypes.INSERTION or variant_call.variant_type == VariantTypes.DELETION:
            alt_tumor_reads = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index('TIR'),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            alt_normal_reads = safely_retrieve_value(
                dict=normal_sample,
                key=format.index('TIR'),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=str
            ).split(',')[0]
            variant_call.alt_tumor_reads = safely_convert_value(
                value=alt_tumor_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
            variant_call.alt_normal_reads = safely_convert_value(
                value=alt_normal_reads,
                type=int,
                default_value=DEFAULT_ATTRIBUTE_VALUE
            )
            variant_call.ref_tumor_reads = variant_call.total_tumor_reads - variant_call.alt_tumor_reads
            variant_call.ref_normal_reads = variant_call.total_normal_reads - variant_call.alt_normal_reads

        # Update total_tumor_reads if it is currently unknown but can be inferred
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.ref_tumor_reads) == int and \
            variant_call.total_tumor_reads is None:
            variant_call.total_tumor_reads = variant_call.alt_tumor_reads + variant_call.ref_tumor_reads

        # Update total_normal_reads if it is currently unknown but can be inferred
        if type(variant_call.alt_normal_reads) == int and \
            type(variant_call.ref_normal_reads) == int and \
            variant_call.total_normal_reads is None:
            variant_call.total_normal_reads = variant_call.alt_normal_reads + variant_call.ref_normal_reads

        # Update alt_tumor_fraction if it is currently unknown
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.total_tumor_reads) == int and \
            variant_call.alt_tumor_fraction is None:
            if variant_call.total_tumor_reads > 0:
                variant_call.alt_tumor_fraction = float(variant_call.alt_tumor_reads) / float(variant_call.total_tumor_reads)

        # Update alt_normal_fraction if it is currently unknown
        if type(variant_call.alt_normal_reads) == int and \
            type(variant_call.total_normal_reads) == int and \
            variant_call.alt_normal_fraction is None:
            if variant_call.total_normal_reads > 0:
                variant_call.alt_normal_fraction = float(variant_call.alt_normal_reads) / float(variant_call.total_normal_reads)

        # Get variant call ID
        variant_call.id = '%s_%s_%s' % \
                          (VariantCallingMethods.STRELKA2_SOMATIC,
                           variant_call.variant_type,
                           curr_idx)
        variant = Variant(id='variant_%i' % curr_idx)
        variant.variant_calls[variant_call.id] = variant_call
        variants_list.variants[variant.id] = variant
        curr_idx += 1

    logger.info('%i rows in the returning VariantsList.' % len(variants_list.variants))
    return variants_list


def parse_strelka2_germline_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str
    ) -> VariantsList:
    """
    Parses a Strelka2-Germline DataFrame and
    returns an instance of the VariantsList class.

    Parameters
    ----------
    df_vcf                  :   DataFrame of rows from a Strelka2 (germline mode) VCF file.
    sequencing_platform     :   Sequencing platform.
    source_id               :   Source ID.

    Returns
    -------
    variants_list           :   An instance of the VariantsList class.
    """
    variants_list = VariantsList()
    tumor_sample_id = df_vcf.columns.values.tolist()[-1]
    curr_idx = 1
    for row in df_vcf.to_dict('records'):
        variant_call = VariantCall()
        variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_id'] = safely_retrieve_value(
            dict=row,
            key='ID',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.nucleic_acid = NucleicAcidTypes.DNA
        variant_call.source_id = source_id
        variant_call.tumor_sample_id = tumor_sample_id
        variant_call.variant_calling_method = VariantCallingMethods.STRELKA2_GERMLINE
        variant_call.sequencing_platform = sequencing_platform
        variant_call.chr_1 = safely_retrieve_value(
            dict=row,
            key='CHROM',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.pos_1 = safely_retrieve_value(
            dict=row,
            key='POS',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=int
        )
        variant_call.chr_2 = safely_retrieve_value(
            dict=row,
            key='CHROM',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.ref = safely_retrieve_value(
            dict=row,
            key='REF',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.alt = safely_retrieve_value(
            dict=row,
            key='ALT',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.filter = safely_retrieve_value(
            dict=row,
            key='FILTER',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.quality_score = safely_retrieve_value(
            dict=row,
            key='QUAL',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=float
        )

        if len(variant_call.ref) == 1 and len(variant_call.alt) == 1:
            variant_call.variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
            variant_call.variant_sequences.append(variant_call.alt)
            variant_call.variant_size = 1
            variant_call.pos_2 = variant_call.pos_1
        elif len(variant_call.ref) == 1 and len(variant_call.alt) > 1:
            variant_call.variant_type = VariantTypes.INSERTION
            variant_call.variant_sequences.append(variant_call.alt[1:])
            variant_call.variant_size = len(variant_call.alt[1:])
            variant_call.pos_2 = variant_call.pos_1
        elif len(variant_call.ref) > 1 and len(variant_call.alt) == 1:
            variant_call.variant_type = VariantTypes.DELETION
            variant_call.variant_sequences.append(variant_call.ref[1:])
            variant_call.variant_size = len(variant_call.ref[1:])
            variant_call.pos_1 = variant_call.pos_1 + 1
            variant_call.pos_2 = variant_call.pos_1 + variant_call.variant_size - 1
        elif len(variant_call.ref) > 1 and len(variant_call.alt) > 1:
            variant_call.variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
            variant_call.variant_sequences.append(variant_call.alt)
            variant_call.variant_size = len(variant_call.alt)
            variant_call.pos_2 = variant_call.pos_1 + len(variant_call.variant_sequences[0]) - 1
        else:
            logger.error('Unknown variant type. REF: %s. ALT: %s' % (variant_call.ref, variant_call.alt))
            exit(1)

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            if '=' in curr_info:
                curr_info_elements = curr_info.split('=')
                variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_' + curr_info_elements[0].lower()] = safely_convert_value(
                    value=curr_info_elements[1],
                    default_value=DEFAULT_ATTRIBUTE_VALUE,
                    type=VariantCallingMethods.AttributeTypes.STRELKA2_GERMLINE[VariantCallingMethods.STRELKA2_GERMLINE + '_' + curr_info_elements[0].lower()]
                )
            else:
                variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_' + curr_info.lower()] = True

        # Extract FORMAT
        format = str(row['FORMAT']).split(':')
        tumor_sample = str(row[tumor_sample_id]).split(':')
        for curr_format in format:
            variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_' + curr_format.lower()] = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(curr_format),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=VariantCallingMethods.AttributeTypes.STRELKA2_GERMLINE[VariantCallingMethods.STRELKA2_GERMLINE + '_' + curr_format.lower()]
            )

        # Update total_tumor_reads
        if VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_dp' in variant_call.tool_attributes.keys():
            variant_call.total_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_dp']

        # Update ref_tumor_reads and alt_tumor_reads
        if VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_ad' in variant_call.tool_attributes.keys():
            variant_call.ref_tumor_reads = safely_convert_value(
                value=variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_ad'].split(',')[0],
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=int
            )
            variant_call.alt_tumor_reads = safely_convert_value(
                value=variant_call.tool_attributes[VariantCallingMethods.STRELKA2_GERMLINE + '_tumor_ad'].split(',')[1],
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=int
            )

        # Update total_tumor_reads if it is currently unknown but can be inferred
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.ref_tumor_reads) == int and \
            variant_call.total_tumor_reads is None:
            variant_call.total_tumor_reads = variant_call.alt_tumor_reads + variant_call.ref_tumor_reads

        # Update alt_tumor_fraction if it is currently unknown
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.total_tumor_reads) == int and \
            variant_call.alt_tumor_fraction is None:
            if variant_call.total_tumor_reads > 0:
                variant_call.alt_tumor_fraction = float(variant_call.alt_tumor_reads) / float(variant_call.total_tumor_reads)

        # Get variant call ID
        variant_call.id = '%s_%s_%s' % \
                          (VariantCallingMethods.STRELKA2_GERMLINE,
                           variant_call.variant_type,
                           curr_idx)
        variant = Variant(id='variant_%i' % curr_idx)
        variant.variant_calls[variant_call.id] = variant_call
        variants_list.variants[variant.id] = variant
        curr_idx += 1

        logger.info('%i rows in the returning VariantsList.' % len(variants_list.variants))
        return variants_list