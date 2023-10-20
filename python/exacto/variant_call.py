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
The purpose of this python3 script is to implement the VariantCall dataclass.
"""


import re
import pandas as pd
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import List, Dict
from functools import total_ordering
from .constants import TranslocationOrientations, VariantTypes
from .variant_annotation import VariantAnnotation


@total_ordering
@dataclass(frozen=True)
class VariantCall:
    id: str
    source_id: str
    sample_id: str
    phase_block_id: str
    clone_id: str
    nucleic_acid: str
    variant_calling_method: str
    sequencing_platform: str
    chromosome_1: str
    position_1: int
    chromosome_2: str
    position_2: int
    reference_allele: str
    alternate_allele: str
    filter: str
    quality_score: float
    precise: bool
    variant_type: str
    variant_subtype: str
    variant_size: int
    total_read_count: int
    reference_allele_read_count: int
    alternate_allele_read_count: int
    alternate_allele_fraction: float
    alternate_allele_read_ids: List[str] = field(default_factory=list)
    variant_sequences: List[str] = field(default_factory=list)
    tool_attributes: OrderedDict = field(default_factory=dict)
    position_1_annotations: List[VariantAnnotation] = field(default_factory=list)
    position_2_annotations: List[VariantAnnotation] = field(default_factory=list)
    tags: List[str] = field(default_factory=list)

    def __lt__(self, other):
        if isinstance(other, VariantCall):
            return (self.chromosome_1,
                    self.position_1,
                    self.chromosome_2,
                    self.position_2) < \
                   (other.chromosome_1,
                    other.position_1,
                    other.chromosome_2,
                    other.position_2)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, VariantCall):
            return (self.chromosome_1,
                    self.position_1,
                    self.chromosome_2,
                    self.position_2) == \
                   (other.chromosome_1,
                    other.position_1,
                    other.chromosome_2,
                    other.position_2)
        return NotImplemented

    def add_position_1_annotation(self, variant_annotation: VariantAnnotation):
        self.position_1_annotations.append(variant_annotation)

    def add_position_2_annotation(self, variant_annotation: VariantAnnotation):
        self.position_2_annotations.append(variant_annotation)

    def get_translocation_orientation(self):
        """
        Returns translocation metadata.

        Returns
        -------
        metadata    :   Tuple (orientation, t_chromosome, t_position, p_chromosome, p_position).
        """
        if self.variant_type == VariantTypes.TRANSLOCATION:
            if re.search("^.*\[.*\[$", self.alternate_allele):                  # t[p[ piece extending to the right of p is joined after t
                orientation = TranslocationOrientations.ORIENTATION_1
                alternate_allele_elements = self.alternate_allele.split('[')
                t = alternate_allele_elements[0]
                p = alternate_allele_elements[1]
            elif re.search("^.*\].*\]$", self.alternate_allele):                # t]p] reverse comp piece extending left of p is joined after t
                orientation = TranslocationOrientations.ORIENTATION_2
                alternate_allele_elements = self.alternate_allele.split(']')
                t = alternate_allele_elements[0]
                p = alternate_allele_elements[1]
            elif re.search("^\].*\].*$", self.alternate_allele):                # ]p]t piece extending to the left of p is joined before t
                orientation = TranslocationOrientations.ORIENTATION_3
                alternate_allele_elements = self.alternate_allele.split(']')
                t = alternate_allele_elements[2]
                p = alternate_allele_elements[1]
            elif re.search("^\[.*\[.*$", self.alternate_allele):                # [p[t  reverse comp piece extending right of p is joined before t
                orientation = TranslocationOrientations.ORIENTATION_4
                alternate_allele_elements = self.alternate_allele.split('[')
                t = alternate_allele_elements[2]
                p = alternate_allele_elements[1]
            else:
                raise Exception('Unknown ALT format to infer translocation orientation type: %s' % self.alternate_allele)

            if p == '%s:%i' % (self.chromosome_1, self.position_1):
                p_chromosome = self.chromosome_1
                p_position = self.position_1
                t_chromosome = self.chromosome_2
                t_position = self.position_2
            elif p == '%s:%i' % (self.chromosome_2, self.position_2):
                p_chromosome = self.chromosome_2
                p_position = self.position_2
                t_chromosome = self.chromosome_1
                t_position = self.position_1
            else:
                raise Exception('Positions for p and t could not be inferred from self.alternate_allele: %s' % self.alternate_allele)
            return orientation, t_chromosome, t_position, p_chromosome, p_position
        else:
            raise Exception('This VariantCall object does not encode a translocation. '
                            'Therefore translocation orientation cannot be inferred.')

    def to_dict(self) -> Dict:
        data = {
            'id': self.id,
            'source_id': self.id,
            'sample_id': self.sample_id,
            'phase_block_id': self.phase_block_id,
            'clone_id': self.clone_id,
            'nucleic_acid': self.nucleic_acid,
            'variant_calling_method': self.variant_calling_method,
            'sequencing_platform': self.sequencing_platform,
            'chromosome_1': self.chromosome_1,
            'position_1': self.position_1,
            'chromosome_2': self.chromosome_2,
            'position_2': self.position_2,
            'reference_allele': self.reference_allele,
            'alternate_allele': self.alternate_allele,
            'filter': self.filter,
            'quality_score': self.quality_score,
            'precise': self.precise,
            'variant_type': self.variant_type,
            'variant_subtype': self.variant_subtype,
            'variant_size': self.variant_size,
            'total_read_count': self.total_read_count,
            'reference_allele_read_count': self.reference_allele_read_count,
            'alternate_allele_read_count': self.alternate_allele_read_count,
            'alternate_allele_fraction': self.alternate_allele_fraction,
            'alternate_allele_read_ids': self.alternate_allele_read_ids,
            'variant_sequences': self.variant_sequences
        }
        tool_attributes_ = {}
        for key, value in self.tool_attributes.items():
            tool_attributes_[key] = str(value)
        data['tool_attributes'] = tool_attributes_
        data['position_1_annotations'] = [attribute.to_dict() for attribute in self.position_1_annotations]
        data['position_2_annotations'] = [attribute.to_dict() for attribute in self.position_2_annotations]
        data['tags'] = self.tags
        return data

    def to_dataframe_row(self) -> Dict:
        data = {
            'variant_call_id': ['' if self.id is None else self.id],
            'source_id': ['' if self.source_id is None else self.source_id],
            'sample_id': ['' if self.sample_id is None else self.sample_id],
            'phase_block_id': ['' if self.phase_block_id is None else self.phase_block_id],
            'clone_id': ['' if self.clone_id is None else self.clone_id],
            'nucleic_acid': ['' if self.nucleic_acid is None else self.nucleic_acid],
            'variant_calling_method': ['' if self.variant_calling_method is None else self.variant_calling_method],
            'sequencing_platform': ['' if self.sequencing_platform is None else self.sequencing_platform],
            'chromosome_1': ['' if self.chromosome_1 is None else self.chromosome_1],
            'position_1': ['' if self.position_1 is None else self.position_1],
            'chromosome_2': ['' if self.chromosome_2 is None else self.chromosome_2],
            'position_2': ['' if self.position_2 is None else self.position_2],
            'reference_allele': ['' if self.reference_allele is None else self.reference_allele],
            'alternate_allele': ['' if self.alternate_allele is None else self.alternate_allele],
            'filter': ['' if self.filter is None else self.filter],
            'quality_score': ['' if self.quality_score is None else self.quality_score],
            'precise': ['' if self.precise is None else self.precise],
            'variant_type': ['' if self.variant_type is None else self.variant_type],
            'variant_subtype': ['' if self.variant_subtype is None else self.variant_subtype],
            'variant_size': ['' if self.variant_size is None else self.variant_size],
            'variant_sequences': [';'.join(self.variant_sequences)],
            'total_read_count': ['' if self.total_read_count is None else self.total_read_count],
            'reference_allele_read_count': ['' if self.reference_allele_read_count is None else self.reference_allele_read_count],
            'alternate_allele_read_count': ['' if self.alternate_allele_read_count is None else self.alternate_allele_read_count],
            'alternate_allele_fraction': ['' if self.alternate_allele_fraction is None else self.alternate_allele_fraction],
            'alternate_allele_read_ids': [';'.join(self.alternate_allele_read_ids)]
        }

        tool_attributes = []
        for key, val in self.tool_attributes.items():
            tool_attributes.append('%s=%s' % (key, val))
        data['tool_attributes'] = [';'.join(tool_attributes)]

        pos_1_annotation_region = []
        pos_1_annotation_source = []
        pos_1_annotation_source_version = []
        pos_1_annotation_gene_id = []
        pos_1_annotation_gene_stable_id = []
        pos_1_annotation_gene_version = []
        pos_1_annotation_gene_name = []
        pos_1_annotation_gene_type = []
        pos_1_annotation_gene_strand = []
        pos_1_annotation_species = []
        for i in self.position_1_annotations:
            pos_1_annotation_region.append(i.region)
            pos_1_annotation_source.append(i.source)
            pos_1_annotation_source_version.append(i.source_version)
            pos_1_annotation_gene_id.append('' if i.gene_id is None else i.gene_id)
            pos_1_annotation_gene_stable_id.append('' if i.gene_stable_id is None else i.gene_stable_id)
            pos_1_annotation_gene_version.append('' if i.gene_version is None else i.gene_version)
            pos_1_annotation_gene_name.append('' if i.gene_name is None else i.gene_name)
            pos_1_annotation_gene_type.append('' if i.gene_type is None else i.gene_type)
            pos_1_annotation_gene_strand.append('' if i.gene_strand is None else i.gene_strand)
            pos_1_annotation_species.append('' if i.species is None else i.species)

        pos_2_annotation_region = []
        pos_2_annotation_source = []
        pos_2_annotation_source_version = []
        pos_2_annotation_gene_id = []
        pos_2_annotation_gene_stable_id = []
        pos_2_annotation_gene_version = []
        pos_2_annotation_gene_name = []
        pos_2_annotation_gene_type = []
        pos_2_annotation_gene_strand = []
        pos_2_annotation_species = []
        for i in self.position_1_annotations:
            pos_2_annotation_region.append(i.region)
            pos_2_annotation_source.append(i.source)
            pos_2_annotation_source_version.append(i.source_version)
            pos_2_annotation_gene_id.append('' if i.gene_id is None else i.gene_id)
            pos_2_annotation_gene_stable_id.append('' if i.gene_stable_id is None else i.gene_stable_id)
            pos_2_annotation_gene_version.append('' if i.gene_version is None else i.gene_version)
            pos_2_annotation_gene_name.append('' if i.gene_name is None else i.gene_name)
            pos_2_annotation_gene_type.append('' if i.gene_type is None else i.gene_type)
            pos_2_annotation_gene_strand.append('' if i.gene_strand is None else i.gene_strand)
            pos_2_annotation_species.append('' if i.species is None else i.species)

        data['position_1_annotation_region'] = [';'.join(pos_1_annotation_region)]
        data['position_1_annotation_source'] = [';'.join(pos_1_annotation_source)]
        data['position_1_annotation_source_version'] = [';'.join(pos_1_annotation_source_version)]
        data['position_1_annotation_gene_id'] = [';'.join(pos_1_annotation_gene_id)]
        data['position_1_annotation_gene_stable_id'] = [';'.join(pos_1_annotation_gene_stable_id)]
        data['position_1_annotation_gene_version'] = [';'.join(pos_1_annotation_gene_version)]
        data['position_1_annotation_gene_name'] = [';'.join(pos_1_annotation_gene_name)]
        data['position_1_annotation_gene_type'] = [';'.join(pos_1_annotation_gene_type)]
        data['position_1_annotation_gene_strand'] = [';'.join(pos_1_annotation_gene_strand)]
        data['position_1_annotation_species'] = [';'.join(pos_1_annotation_species)]

        data['position_2_annotation_region'] = [';'.join(pos_2_annotation_region)]
        data['position_2_annotation_source'] = [';'.join(pos_2_annotation_source)]
        data['position_2_annotation_source_version'] = [';'.join(pos_2_annotation_source_version)]
        data['position_2_annotation_gene_id'] = [';'.join(pos_2_annotation_gene_id)]
        data['position_2_annotation_gene_stable_id'] = [';'.join(pos_2_annotation_gene_stable_id)]
        data['position_2_annotation_gene_version'] = [';'.join(pos_2_annotation_gene_version)]
        data['position_2_annotation_gene_name'] = [';'.join(pos_2_annotation_gene_name)]
        data['position_2_annotation_gene_type'] = [';'.join(pos_2_annotation_gene_type)]
        data['position_2_annotation_gene_strand'] = [';'.join(pos_2_annotation_gene_strand)]
        data['position_2_annotation_species'] = [';'.join(pos_2_annotation_species)]

        data['tags'] = [';'.join(self.tags)]
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())
