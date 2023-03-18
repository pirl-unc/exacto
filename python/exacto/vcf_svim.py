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
related to handling svim VCF files.
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


def parse_svim_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str
    ) -> VariantsList:
    """
    Parses a SVIM DataFrame and returns an instance of the VariantsList class.

    Parameters
    ----------
    df_vcf                  :   DataFrame of rows from a SVIM VCF file.
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
        variant_call.tool_attributes[VariantCallingMethods.SVIM + '_id'] = safely_retrieve_value(
            dict=row,
            key='ID',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=str
        )
        variant_call.nucleic_acid = NucleicAcidTypes.DNA
        variant_call.source_id = source_id
        variant_call.tumor_sample_id = tumor_sample_id
        variant_call.variant_calling_method = VariantCallingMethods.SVIM
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

        # Extract INFO
        info = str(row['INFO']).split(';')
        for curr_info in info:
            if '=' in curr_info:
                curr_info_elements = curr_info.split('=')
                variant_call.tool_attributes[VariantCallingMethods.SVIM + '_' + curr_info_elements[0].lower()] = safely_convert_value(
                    value=curr_info_elements[1],
                    default_value=DEFAULT_ATTRIBUTE_VALUE,
                    type=VariantCallingMethods.AttributeTypes.SVIM[VariantCallingMethods.SVIM + '_' + curr_info_elements[0].lower()]
                )
            else:
                variant_call.tool_attributes[VariantCallingMethods.PBSV + '_' + curr_info.lower()] = True

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        tumor_sample = str(row[tumor_sample_id]).split(':')
        for curr_format in format:
            variant_call.tool_attributes[VariantCallingMethods.SVIM + '_' + curr_format.lower()] = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(curr_format),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=VariantCallingMethods.AttributeTypes.SVIM[VariantCallingMethods.SVIM + '_' + curr_format.lower()]
            )

        # Update variant_type
        if VariantCallingMethods.SVIM + '_svtype' in variant_call.tool_attributes.keys():
            variant_call.variant_type = variant_call.tool_attributes[VariantCallingMethods.SVIM + '_svtype']

        if variant_call.variant_type == 'DUP:TANDEM':
            variant_call.variant_type = VariantTypes.DUPLICATION
            variant_call.variant_subtype = VariantTypes.DuplicationSubtypes.TANDEM_DUPLICATION

        # Update pos_2
        if VariantCallingMethods.SVIM + '_end' in variant_call.tool_attributes.keys():
            variant_call.pos_2 = variant_call.tool_attributes[VariantCallingMethods.SVIM + '_end']

        # Update ref_tumor_reads and alt_tumor_reads
        if VariantCallingMethods.SVIM + '_ad' in variant_call.tool_attributes.keys():
            variant_call.ref_tumor_reads = safely_convert_value(
                value=variant_call.tool_attributes[VariantCallingMethods.SVIM + '_ad'].split(',')[0],
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=int
            )
            variant_call.alt_tumor_reads = safely_convert_value(
                value=variant_call.tool_attributes[VariantCallingMethods.SVIM + '_ad'].split(',')[1],
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=int
            )

        # Update alt_tumor_reads
        if VariantCallingMethods.SVIM + '_support' in variant_call.tool_attributes.keys():
            variant_call.alt_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.SVIM + '_support']

        # Update variant_size
        if VariantCallingMethods.SVIM + '_svlen' in variant_call.tool_attributes.keys():
            variant_call.variant_size = abs(variant_call.tool_attributes[VariantCallingMethods.SVIM + '_svlen'])

        # Update variant_sequences
        if VariantCallingMethods.SVIM + '_seqs' in variant_call.tool_attributes.keys():
            variant_call.variant_sequences = variant_call.tool_attributes[VariantCallingMethods.SVIM + '_seqs'].split(',')

        # Update alt_tumor_read_ids
        if VariantCallingMethods.SVIM + '_reads' in variant_call.tool_attributes.keys():
            variant_call.alt_tumor_read_ids = variant_call.tool_attributes[VariantCallingMethods.SVIM + '_reads'].split(',')

        # Update chr_2 for 'BND'
        if variant_call.variant_type == VariantTypes.BREAKPOINT or variant_call.variant_type == VariantTypes.TRANSLOCATION:
            alt_val = str(row['ALT']).split(":")[0]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            variant_call.chr_2 = str(alt_val)

        # Update pos_2 for 'BND'
        if variant_call.variant_type == VariantTypes.BREAKPOINT or variant_call.variant_type == VariantTypes.TRANSLOCATION:
            alt_val = str(row['ALT']).split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            variant_call.pos_2 = int(alt_val)

        # Update variant_size 'BND'
        if (variant_call.variant_type == VariantTypes.BREAKPOINT) or (variant_call.variant_type == VariantTypes.TRANSLOCATION) and \
                (variant_call.chr_1 == variant_call.chr_2):
            variant_call.variant_size = abs(variant_call.pos_2 - variant_call.pos_1)

        # Update total_tumor_reads if it is currently unknown but can be inferred
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.ref_tumor_reads) == int and \
            variant_call.total_tumor_reads is None:
            variant_call.total_tumor_reads = variant_call.alt_tumor_reads + variant_call.ref_tumor_reads

        # Update alt_tumor_fraction if it is currently unknown but can be inferred
        if type(variant_call.alt_tumor_reads) == int and \
            type(variant_call.total_tumor_reads) == int and \
            variant_call.alt_tumor_fraction is None:
            if variant_call.total_tumor_reads > 0:
                variant_call.alt_tumor_fraction = float(variant_call.alt_tumor_reads) / float(variant_call.total_tumor_reads)

        # Get variant call ID
        variant_call.id = '%s_%s_%s' % \
                          (VariantCallingMethods.SVIM,
                           variant_call.variant_type,
                           curr_idx)
        variant = Variant(id='variant_%i' % curr_idx)
        variant.variant_calls.append(variant_call)
        variants_list.variants.append(variant)
        curr_idx += 1

    logger.info('%i rows in the returning VariantsList.' % len(variants_list.variants))
    return variants_list
