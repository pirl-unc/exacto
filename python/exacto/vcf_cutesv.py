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
related to handling cuteSV VCF files.
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


def parse_cutesv_callset(
        df_vcf: pd.DataFrame,
        sequencing_platform: str,
        source_id: str
    ) -> VariantsList:
    """
    Parses a cuteSV DataFrame and
    returns an instance of the VariantsList class.

    Parameters
    ----------
    df_vcf                  :   DataFrame of rows from a cuteSV VCF file.
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
        variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_id'] = safely_retrieve_value(
            dict=row,
            key='ID',
            default_value=DEFAULT_ATTRIBUTE_VALUE,
            type=VariantCallingMethods.AttributeTypes.CUTESV[VariantCallingMethods.CUTESV + '_id']
        )
        variant_call.nucleic_acid = NucleicAcidTypes.DNA
        variant_call.source_id = source_id
        variant_call.tumor_sample_id = tumor_sample_id
        variant_call.variant_calling_method = VariantCallingMethods.CUTESV
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
                variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_' + curr_info_elements[0].lower()] = safely_convert_value(
                    value=curr_info_elements[1],
                    default_value=DEFAULT_ATTRIBUTE_VALUE,
                    type=VariantCallingMethods.AttributeTypes.CUTESV[VariantCallingMethods.CUTESV + '_' + curr_info_elements[0].lower()]
                )
            else:
                if curr_info == 'PRECISE':
                    variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_precise'] = True
                if curr_info == 'IMPRECISE':
                    variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_precise'] = False

        # Extract FORMAT (sample)
        format = str(row['FORMAT']).split(':')
        tumor_sample = str(row[tumor_sample_id]).split(':')
        for curr_format in format:
            variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_' + curr_format.lower()] = safely_retrieve_value(
                dict=tumor_sample,
                key=format.index(curr_format),
                default_value=DEFAULT_ATTRIBUTE_VALUE,
                type=VariantCallingMethods.AttributeTypes.CUTESV[VariantCallingMethods.CUTESV + '_' + curr_format.lower()]
            )

        # Update variant_type
        if VariantCallingMethods.CUTESV + '_svtype' in variant_call.tool_attributes.keys():
            variant_call.variant_type = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_svtype']

        # Update variant_size
        if VariantCallingMethods.CUTESV + '_svlen' in variant_call.tool_attributes.keys():
            variant_call.variant_size = abs(variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_svlen'])

        # Update precise
        if VariantCallingMethods.CUTESV + '_precise' in variant_call.tool_attributes.keys():
            variant_call.precise = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_precise']

        # Update chr_2
        if VariantCallingMethods.CUTESV + '_chr2' in variant_call.tool_attributes.keys():
            variant_call.chr_2 = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_chr2']

        # Update pos_2
        if VariantCallingMethods.CUTESV + '_end' in variant_call.tool_attributes.keys():
            variant_call.pos_2 = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_end']

        # Update ref_tumor_reads
        if VariantCallingMethods.CUTESV + '_dr' in variant_call.tool_attributes.keys():
            variant_call.ref_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_dr']

        # Update alt_tumor_reads
        if VariantCallingMethods.CUTESV + '_re' in variant_call.tool_attributes.keys():
            variant_call.alt_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_re']
        if variant_call.alt_tumor_reads is None:
            variant_call.alt_tumor_reads = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_dv']

        # Update alt_tumor_read_ids
        if VariantCallingMethods.CUTESV + '_rnames' in variant_call.tool_attributes.keys():
            variant_call.alt_tumor_read_ids = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_rnames'].split(',')

        # Update alt_tumor_fraction
        if VariantCallingMethods.CUTESV + '_af' in variant_call.tool_attributes.keys():
            variant_call.alt_tumor_fraction = variant_call.tool_attributes[VariantCallingMethods.CUTESV + '_af']

        # Update chr_2 for 'BND' or 'TRA'
        if variant_call.variant_type == VariantTypes.BREAKPOINT or variant_call.variant_type == VariantTypes.TRANSLOCATION:
            chr_2 = variant_call.alt.split(":")[0]
            chr_2 = chr_2.replace("[", "")
            chr_2 = chr_2.replace("]", "")
            variant_call.chr_2 = chr_2.replace("N", "")

        # Update pos_2 for 'BND' or 'TRA'
        if variant_call.variant_type == VariantTypes.BREAKPOINT or variant_call.variant_type == VariantTypes.TRANSLOCATION:
            pos_2 = variant_call.alt.split(":")[1]
            pos_2 = pos_2.replace("[", "")
            pos_2 = pos_2.replace("]", "")
            pos_2 = pos_2.replace("N", "")
            variant_call.pos_2 = int(pos_2)

        # Update variant_size for 'BND' or 'TRA'
        if (variant_call.variant_type == VariantTypes.BREAKPOINT or variant_call.variant_type == VariantTypes.TRANSLOCATION) and \
            (variant_call.chr_1 == variant_call.chr_2):
            variant_call.variant_size = abs(variant_call.pos_2 - variant_call.pos_1)

        # Update variant_sequence
        if variant_call.variant_type == VariantTypes.INSERTION:
            variant_call.variant_sequences.append(variant_call.alt[1:])

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

        # Update variant_call ID
        variant_call.id = '%s_%s_%s' % \
                          (VariantCallingMethods.CUTESV,
                           variant_call.variant_type,
                           curr_idx)
        variant = Variant(id='variant_%i' % curr_idx)
        variant.variant_calls.append(variant_call)
        variants_list.variants.append(variant)
        curr_idx += 1

    logger.info('%i rows in the returning VariantsList.' % len(variants_list.variants))
    return variants_list
