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
The purpose of this python3 script is to implement the VCF dataclass.
"""


import gzip
import pandas as pd
from collections import OrderedDict
from .constants import *
from .logging import get_logger
from .utilities import retrieve_with_default, get_typed_value
from .variant_call import VariantCall
from .variant import Variant
from .variants_list import VariantsList


logger = get_logger(__name__)


class Vcf:

    @staticmethod
    def read_vcf_file(vcf_file: str) -> pd.DataFrame:
        """
        Reads a VCF file and returns a DataFrame.

        Parameters
        ----------
        vcf_file    :   VCF file.

        Returns
        -------
        df_vcf      :   DataFrame with the following columns:
                        'CHROM'
                        'POS'
                        'ID'
                        'REF'
                        'ALT'
                        'QUAL'
                        'FILTER'
                        'INFO'
                        'FORMAT'
                        Sample 1
                        Sample 2
                        ...
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
            return pd.read_csv(vcf_file,
                               compression='gzip',
                               comment='#',
                               delim_whitespace=True,
                               header=None,
                               low_memory=True,
                               memory_map=False,
                               names=vcf_names)
        else:
            return pd.read_csv(vcf_file,
                               comment='#',
                               delim_whitespace=True,
                               header=None,
                               low_memory=True,
                               memory_map=False,
                               names=vcf_names)

    @staticmethod
    def parse_cutesv_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a cuteSV DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a cuteSV VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.CUTESV[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        else:
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.CUTESV[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_size
                # chromosome_2
                # position_2
                # reference_allele_read_count
                # alternate_allele_read_count
                # alternate_allele_read_ids
                # alternate_allele_fraction
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])
                if 'precise' in tool_attributes.keys():
                    precise = tool_attributes['precise']
                if 'chr2' in tool_attributes.keys():
                    chromosome_2 = tool_attributes['chr2']
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'dr' in tool_attributes.keys():
                    reference_allele_read_count = tool_attributes['dr']
                if 're' in tool_attributes.keys():
                    alternate_allele_read_count = tool_attributes['re']
                if alternate_allele_read_count is None:
                    alternate_allele_read_count = tool_attributes['dv']
                if 'rnames' in tool_attributes.keys():
                    alternate_allele_read_ids = tool_attributes['rnames'].split(',')
                if 'af' in tool_attributes.keys():
                    alternate_allele_fraction = tool_attributes['af']

                # Update chromosome_2 for 'BND' or 'TRA'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    chromosome_2 = alternate_allele.split(":")[0]
                    chromosome_2 = chromosome_2.replace("[", "")
                    chromosome_2 = chromosome_2.replace("]", "")
                    chromosome_2 = chromosome_2.replace("N", "")

                # Update position_2 for 'BND' or 'TRA'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    position_2 = alternate_allele.split(":")[1]
                    position_2 = position_2.replace("[", "")
                    position_2 = position_2.replace("]", "")
                    position_2 = position_2.replace("N", "")
                    position_2 = int(position_2)

                # Update variant_size for 'BND' or 'TRA'
                if (variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_2 - position_1)

                # Update variant_sequence for 'INS'
                if variant_type == VariantTypes.INSERTION:
                    variant_sequences.append(str(alternate_allele[1:]))

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.CUTESV,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.CUTESV,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    precise=precise,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_deepvariant_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a DeepVariant DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a DeepVariant VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = True
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_subtype = ''
                variant_sequences = []
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Update the following variables:
                # variant_type
                # variant_sequences
                # variant_size
                # position_1
                # position_2
                if len(reference_allele) == 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = 1
                    position_2 = position_1
                elif len(reference_allele) == 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.INSERTION
                    variant_sequences.append(str(alternate_allele[1:]))
                    variant_size = len(alternate_allele[1:])
                    position_2 = position_1
                elif len(reference_allele) > 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.DELETION
                    variant_sequences.append(str(reference_allele[1:]))
                    variant_size = len(reference_allele[1:])
                    position_1 = position_1 + 1
                    position_2 = position_1 + variant_size - 1
                elif len(reference_allele) > 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = len(alternate_allele)
                    position_2 = position_1 + len(variant_sequences[0]) - 1
                else:
                    raise Exception(
                        'Unknown variant type. REF: %s. ALT: %s' %
                        (reference_allele, alternate_allele)
                    )

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.DEEPVARIANT[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        elif curr_info != '.':
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.DEEPVARIANT[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # total_read_count
                # reference_allele_read_count
                # alternate_allele_read_count
                # alternate_allele_fraction
                if 'dp' in tool_attributes.keys():
                    total_read_count = tool_attributes['dp']
                if 'ad' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[1], default_value=-1, type=int)
                if 'vaf' in tool_attributes.keys():
                    if tool_attributes['vaf'] == '':
                        alternate_allele_fraction = -1.0
                    else:
                        alternate_allele_fraction = float(tool_attributes['vaf'])

                # Update total_read_count if it is currently unknown but can be inferred
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count

                # Update alternate_allele_fraction if it is currently unknown but can be inferred
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.DEEPVARIANT,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.DEEPVARIANT,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_gatk4_mutect2_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a GATK4-Mutect2 DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a GATK4-Mutect2 VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = True
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.GATK4_MUTECT2[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.GATK4_MUTECT2[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_sequences
                # variant_size
                # position_1
                # position_2
                if len(reference_allele) == 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = 1
                    position_2 = position_1
                elif len(reference_allele) == 1 and len(alternate_allele) > 1:
                    if ',' in alternate_allele:
                        variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                        variant_sequences.append(str(alternate_allele))
                        position_2 = position_1
                    else:
                        variant_type = VariantTypes.INSERTION
                        variant_sequences.append(str(alternate_allele[1:]))
                        variant_size = len(alternate_allele[1:])
                        position_2 = position_1
                elif len(reference_allele) > 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.DELETION
                    variant_sequences.append(str(reference_allele[1:]))
                    variant_size = len(reference_allele[1:])
                    position_1 = position_1 + 1
                    position_2 = position_1 + variant_size - 1
                elif len(reference_allele) > 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = len(alternate_allele)
                    position_2 = position_1 + len(variant_sequences[0]) - 1
                else:
                    raise Exception(
                        'Unknown variant type. REF: %s. ALT: %s' %
                        (reference_allele, alternate_allele)
                    )

                # Update the following variables:
                # reference_allele_read_count
                # alternate_allele_read_count
                # total_read_count
                if 'ad' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[1], default_value=-1, type=int)
                if 'dp' in tool_attributes.keys():
                    total_read_count = get_typed_value(value=tool_attributes['dp'], default_value=-1, type=int)

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.GATK4_MUTECT2,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.GATK4_MUTECT2,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_pbsv_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a pbsv DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a pbsv VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        included_mate_ids = set()
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.PBSV[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        else:
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.PBSV[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # precise
                # position_2
                # total_read_count
                # reference_allele_read_count
                # alternate_allele_read_count
                # variant_size
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if 'precise' in tool_attributes.keys():
                    precise = tool_attributes['precise']
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'dp' in tool_attributes.keys():
                    total_read_count = tool_attributes['dp']
                if 'ad' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])

                # Update chromosome_2 for 'BND'
                if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                    curr_id = tool_attributes['id'].split("-")[1]
                    chromosome_2 = str(curr_id.split(":")[0])

                # Update position_2 for 'BND'
                if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                    curr_id = tool_attributes['id'].split("-")[1]
                    position_2 = int(curr_id.split(":")[1])

                # Update variant_size for 'BND'
                if (variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_2 - position_1)

                # Update insertion sequence fdor 'INS'
                if variant_type == VariantTypes.INSERTION:
                    variant_sequences.append(str(alternate_allele))

                # Update the following variables if they are currently unknown but can be inferred:
                # total_tumor_reads
                # alt_tumor_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Check if variant_call ID has been included
                if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                    if tool_attributes['id'] in included_mate_ids:
                        continue
                    included_mate_ids.add(tool_attributes['id'])

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.PBSV,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.PBSV,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    precise=precise,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_sniffles2_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a Sniffles2 DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a Sniffles2 VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.SNIFFLES2[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        else:
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.SNIFFLES2[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # precise
                # variant_size
                # variant_type
                # chromosome_2
                # position_2
                # reference_allele_read_count
                # alternate_allele_read_count
                # alternate_allele_fraction
                # alternate_allele_read_ids
                if 'precise' in tool_attributes.keys():
                    precise = tool_attributes['precise']
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if 'chr2' in tool_attributes.keys():
                    chromosome_2 = tool_attributes['chr2']
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'dr' in tool_attributes.keys():
                    reference_allele_read_count = tool_attributes['dr']
                if 'dv' in tool_attributes.keys():
                    alternate_allele_read_count = tool_attributes['dv']
                if 'af' in tool_attributes.keys():
                    alternate_allele_fraction = tool_attributes['af']
                if 'rnames' in tool_attributes.keys():
                    alternate_allele_read_ids = tool_attributes['rnames'].split(',')

                # Update position_2 for 'BND'
                if variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION:
                    alt_val = str(row['ALT']).split(":")[1]
                    alt_val = alt_val.replace("[", "")
                    alt_val = alt_val.replace("]", "")
                    alt_val = alt_val.replace("N", "")
                    position_2 = int(alt_val)

                # Update variant_size for 'BND'
                if (variant_type == VariantTypes.BREAKPOINT or variant_type == VariantTypes.TRANSLOCATION) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_1 - position_2) + 1

                # Update variant_sequence for 'INS'
                if variant_type == VariantTypes.INSERTION:
                    variant_sequences.append(str(alternate_allele))

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.SNIFFLES2,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.SNIFFLES2,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    precise=precise,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_strelka2_somatic_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str,
            case_id: str,
            control_id: str
    ) -> VariantsList:
        """
        Parses a Strelka2 DataFrame and returns VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a Strelka2 (germline mode) VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.
        case_id                 :   Case ID.
        control_id              :   Control ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = True
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.STRELKA2_SOMATIC[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_size
                # variant_sequences
                # position_1
                # position_2
                if len(reference_allele) == 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = 1
                    position_2 = position_1
                elif len(reference_allele) == 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.INSERTION
                    variant_sequences.append(str(alternate_allele[1:]))
                    variant_size = len(alternate_allele[1:])
                    position_2 = position_1
                elif len(reference_allele) > 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.DELETION
                    variant_sequences.append(str(reference_allele[1:]))
                    variant_size = len(reference_allele[1:])
                    position_1 = position_1 + 1
                    position_2 = position_1 + variant_size - 1
                elif len(reference_allele) > 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = len(alternate_allele)
                    position_2 = position_1 + len(variant_sequences[0]) - 1
                else:
                    raise Exception(
                        'Unknown variant type. REF: %s. ALT: %s' %
                        (reference_allele, alternate_allele)
                    )

                # Update the following variables:
                # total_read_count
                # reference_allele_read_count
                # alternate_allele_read_count
                if 'dp' in tool_attributes.keys():
                    total_read_count = tool_attributes['dp']
                if 'ad' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[1], default_value=-1, type=int)

                # Update the following variable if they are currently unknown but can be inferred
                # total_read_count
                # alternate_allele_fraction
                if type(reference_allele_read_count) == int and \
                    type(alternate_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.STRELKA2_SOMATIC,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)

                # Replace sample_id
                if sample_id == 'TUMOR':
                    sample_id = case_id
                if sample_id == 'NORMAL':
                    sample_id = control_id

                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.STRELKA2_SOMATIC,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_svim_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a SVIM DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a SVIM VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.SVIM[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.SVIM[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_subtype
                # position_2
                # reference_allele_read_count
                # alternate_allele_read_count
                # variant_size
                # variant_sequences
                # alternate_allele_read_ids
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if variant_type == 'DUP:TANDEM':
                    variant_type = VariantTypes.DUPLICATION
                    variant_subtype = VariantTypes.DuplicationSubtypes.TANDEM_DUPLICATION
                if variant_type == 'DUP:INT':
                    variant_type = VariantTypes.DUPLICATION
                    variant_subtype = VariantTypes.DuplicationSubtypes.INTERSPERSED_DUPLICATION
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'ad' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[0], default_value=-1, type=int)
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['ad'].split(',')[1], default_value=-1, type=int)
                if (alternate_allele_read_count is None or alternate_allele_read_count == '.' or alternate_allele_read_count == -1) and 'support' in tool_attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['support'], default_value=-1, type=int)
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])
                if 'seqs' in tool_attributes.keys():
                    variant_sequences = [str(sequence) for sequence in tool_attributes['seqs'].split(',')]
                if 'reads' in tool_attributes.keys():
                    alternate_allele_read_ids = tool_attributes['reads'].split(',')

                # Update chromosome_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    alt_val = alternate_allele.split(":")[0]
                    alt_val = alt_val.replace("[", "")
                    alt_val = alt_val.replace("]", "")
                    alt_val = alt_val.replace("N", "")
                    chromosome_2 = str(alt_val)

                # Update position_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    alt_val = alternate_allele.split(":")[1]
                    alt_val = alt_val.replace("[", "")
                    alt_val = alt_val.replace("]", "")
                    alt_val = alt_val.replace("N", "")
                    position_2 = int(alt_val)

                # Update variant_size 'BND'
                if (variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_2 - position_1)

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.SVIM,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.SVIM,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_delly2_somatic_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a Delly2 somatic DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a Delly2 somatic VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.DELLY2_SOMATIC[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        else:
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.DELLY2_SOMATIC[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_subtype
                # position_2
                # reference_allele_read_count
                # alternate_allele_read_count
                # variant_size
                # variant_sequences
                # alternate_allele_read_ids
                if 'precise' in tool_attributes.keys():
                    precise = True
                if 'imprecise' in tool_attributes.keys():
                    precise = False
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'dr' in tool_attributes.keys():
                    reference_allele_read_count = get_typed_value(value=tool_attributes['dr'], default_value=-1, type=int)
                if 'dv' in tool_attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['dv'], default_value=-1, type=int)
                if alternate_allele_read_count is None and 'pe' in tool_attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['pe'], default_value=-1, type=int)
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])
                if 'consensus' in tool_attributes.keys():
                    variant_sequences = [str(tool_attributes['consensus'])]

                # Update chromosome_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    chromosome_2 = str(tool_attributes['chr2'])

                # Update position_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    position_2 = int(tool_attributes['pos2'])

                # Update variant_size 'BND'
                if (variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_2 - position_1)

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.DELLY2_SOMATIC,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.DELLY2_SOMATIC,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_lumpy_somatic_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a Lumpy somatic DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a Lumpy somatic VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.LUMPY_SOMATIC[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        if curr_info == 'PRECISE':
                            tool_attributes['precise'] = True
                        elif curr_info == 'IMPRECISE':
                            tool_attributes['precise'] = False
                        else:
                            tool_attributes[curr_info.lower()] = True

                # Extract FORMAT
                format = str(row['FORMAT']).split(':')
                curr_sample = str(row[sample_id]).split(':')
                for curr_format in format:
                    curr_key = curr_format.lower()
                    curr_type = VariantCallingMethods.AttributeTypes.LUMPY_SOMATIC[curr_key]
                    tool_attributes[curr_key] = retrieve_with_default(dict=curr_sample, key=format.index(curr_format), default_value='', type=curr_type)

                # Update the following variables:
                # variant_type
                # variant_subtype
                # position_2
                # reference_allele_read_count
                # alternate_allele_read_count
                # variant_size
                # variant_sequences
                # alternate_allele_read_ids
                if 'precise' in tool_attributes.keys():
                    precise = True
                if 'imprecise' in tool_attributes.keys():
                    precise = False
                if 'svtype' in tool_attributes.keys():
                    variant_type = tool_attributes['svtype']
                if 'end' in tool_attributes.keys():
                    position_2 = tool_attributes['end']
                if 'su' in tool_attributes.keys():
                    alternate_allele_read_count = get_typed_value(value=tool_attributes['su'], default_value=-1, type=int)
                if 'svlen' in tool_attributes.keys():
                    variant_size = abs(tool_attributes['svlen'])

                # Update chromosome_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    alt_val = alternate_allele.split(":")[0]
                    alt_val = alt_val.replace("[", "")
                    alt_val = alt_val.replace("]", "")
                    alt_val = alt_val.replace("N", "")
                    chromosome_2 = str(alt_val)

                # Update position_2 for 'BND'
                if variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]:
                    alt_val = alternate_allele.split(":")[1]
                    alt_val = alt_val.replace("[", "")
                    alt_val = alt_val.replace("]", "")
                    alt_val = alt_val.replace("N", "")
                    position_2 = int(alt_val)

                # Update variant_size 'BND'
                if (variant_type in [VariantTypes.BREAKPOINT, VariantTypes.TRANSLOCATION]) and \
                    (chromosome_1 == chromosome_2):
                    variant_size = abs(position_2 - position_1)

                # Update the following variables if they are currently unknown but can be inferred:
                # total_read_count
                # alternate_allele_fraction
                if type(alternate_allele_read_count) == int and \
                    type(reference_allele_read_count) == int and \
                    total_read_count is None:
                    total_read_count = alternate_allele_read_count + reference_allele_read_count
                if type(alternate_allele_read_count) == int and \
                    type(total_read_count) == int and \
                    alternate_allele_fraction is None:
                    if total_read_count > 0:
                        alternate_allele_fraction = float(alternate_allele_read_count) / float(total_read_count)

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.LUMPY_SOMATIC,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.LUMPY_SOMATIC,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list

    @staticmethod
    def parse_dbsnp_callset(
            df_vcf: pd.DataFrame,
            sequencing_platform: str,
            source_id: str
    ) -> VariantsList:
        """
        Parses a dbSNP DataFrame and returns a VariantsList object.

        Parameters
        ----------
        df_vcf                  :   DataFrame of rows from a dbSNP VCF file.
        sequencing_platform     :   Sequencing platform.
        source_id               :   Source ID.

        Returns
        -------
        variants_list           :   VariantsList object.
        """
        variants_list = VariantsList()
        sample_ids = df_vcf.columns.values.tolist()[9:]
        curr_variant_call_idx = 1
        curr_variant_idx = 1
        for row in df_vcf.to_dict('records'):
            for sample_id in sample_ids:
                # Initialize values
                phase_block_id = ''
                clone_id = ''
                chromosome_1 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                chromosome_2 = retrieve_with_default(dict=row, key='CHROM', default_value='', type=str)
                position_1 = retrieve_with_default(dict=row, key='POS', default_value=-1, type=int)
                position_2 = -1
                reference_allele = retrieve_with_default(dict=row, key='REF', default_value='', type=str)
                alternate_allele = retrieve_with_default(dict=row, key='ALT', default_value='', type=str)
                filter = retrieve_with_default(dict=row, key='FILTER', default_value='', type=str)
                quality_score = retrieve_with_default(dict=row, key='QUAL', default_value=-1.0, type=float)
                precise = False
                total_read_count = -1
                reference_allele_read_count = -1
                alternate_allele_read_count = -1
                alternate_allele_fraction = -1.0
                variant_type = ''
                variant_subtype = ''
                variant_sequences = []
                variant_size = -1
                alternate_allele_read_ids = []
                tool_attributes = OrderedDict()
                tool_attributes['id'] = retrieve_with_default(dict=row, key='ID', default_value='', type=str)

                # Extract INFO
                info = str(row['INFO']).split(';')
                for curr_info in info:
                    if '=' in curr_info:
                        curr_info_elements = curr_info.split('=')
                        curr_key = curr_info_elements[0].lower()
                        curr_type = VariantCallingMethods.AttributeTypes.CUTESV[curr_key]
                        tool_attributes[curr_key] = get_typed_value(value=curr_info_elements[1], default_value='', type=curr_type)
                    else:
                        tool_attributes[curr_info.lower()] = True

                # Update the following variables:
                # variant_type
                # variant_sequences
                # variant_size
                # position_1
                # position_2
                if len(reference_allele) == 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.SINGLE_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = 1
                    position_2 = position_1
                elif len(reference_allele) == 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.INSERTION
                    variant_sequences.append(str(alternate_allele[1:]))
                    variant_size = len(alternate_allele[1:])
                    position_2 = position_1
                elif len(reference_allele) > 1 and len(alternate_allele) == 1:
                    variant_type = VariantTypes.DELETION
                    variant_sequences.append(str(reference_allele[1:]))
                    variant_size = len(reference_allele[1:])
                    position_1 = position_1 + 1
                    position_2 = position_1 + variant_size - 1
                elif len(reference_allele) > 1 and len(alternate_allele) > 1:
                    variant_type = VariantTypes.MULTI_NUCLEOTIDE_VARIANT
                    variant_sequences.append(str(alternate_allele))
                    variant_size = len(alternate_allele)
                    position_2 = position_1 + variant_sequences[0].length - 1
                else:
                    raise Exception(
                        'Unknown variant type. REF: %s. ALT: %s' %
                        (reference_allele, alternate_allele)
                    )

                # Append variant call to variants list
                variant_call_id = '%s_%s_%s_%s_%s_%s_%i_%s:%i_%s:%i' % (
                    source_id,
                    sample_id,
                    NucleicAcidTypes.DNA,
                    sequencing_platform,
                    VariantCallingMethods.DBSNP,
                    variant_type,
                    curr_variant_call_idx,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2
                )
                variant_id = str(curr_variant_idx)
                variant = Variant(id=variant_id)
                variant_call = VariantCall(
                    id=variant_call_id,
                    source_id=source_id,
                    sample_id=sample_id,
                    phase_block_id=phase_block_id,
                    clone_id=clone_id,
                    nucleic_acid=NucleicAcidTypes.DNA,
                    variant_calling_method=VariantCallingMethods.DBSNP,
                    sequencing_platform=sequencing_platform,
                    chromosome_1=chromosome_1,
                    position_1=position_1,
                    chromosome_2=chromosome_2,
                    position_2=position_2,
                    reference_allele=reference_allele,
                    alternate_allele=alternate_allele,
                    filter=filter,
                    quality_score=quality_score,
                    precise=precise,
                    variant_type=variant_type,
                    variant_subtype=variant_subtype,
                    variant_size=variant_size,
                    variant_sequences=variant_sequences,
                    total_read_count=total_read_count,
                    reference_allele_read_count=reference_allele_read_count,
                    alternate_allele_read_count=alternate_allele_read_count,
                    alternate_allele_fraction=alternate_allele_fraction,
                    alternate_allele_read_ids=alternate_allele_read_ids,
                    tool_attributes=tool_attributes
                )
                variant.add_variant_call(variant_call=variant_call)
                variants_list.add_variant(variant=variant)
                curr_variant_call_idx += 1
            curr_variant_idx += 1

        logger.info('%i variants and %i variant calls in the returning VariantsList.' %
                    (len(variants_list.variant_ids), len(variants_list.variant_call_ids)))
        return variants_list
