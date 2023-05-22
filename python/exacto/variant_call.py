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
The purpose of this python3 script is to implement the VariantCall class.
"""


import re
import pandas as pd
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import List, Dict
from functools import total_ordering
from .constants import TranslocationOrientations, VariantTypes
from .nucleotide_sequence import NucleotideSequence
from .variant_annotation import VariantAnnotation


@total_ordering
@dataclass(frozen=True)
class VariantCall:
    id: str
    source_id: str = None
    sample_id: str = None
    phase_block_id: str = None
    clone_set_id: str = None
    nucleic_acid: str = None
    variant_calling_method: str = None
    sequencing_platform: str = None
    chromosome_1: str = None
    position_1: int = None
    chromosome_2: str = None
    position_2: int = None
    reference_allele: str = None
    alternate_allele: str = None
    filter: str = None
    quality_score: float = None
    precise: bool = None
    variant_type: str = None
    variant_subtype: str = None
    variant_size: int = None
    variant_sequences: List[NucleotideSequence] = field(default_factory=list, repr=False)
    total_read_count: int = None
    reference_allele_read_count: int = None
    alternate_allele_read_count: int = None
    alternate_allele_fraction: float = None
    alternate_allele_read_ids: List[str] = field(default_factory=list)
    tool_attributes: OrderedDict = field(default_factory=dict, repr=False, compare=False)
    position_1_annotations: List[VariantAnnotation] = field(default_factory=list, repr=False, compare=False)
    position_2_annotations: List[VariantAnnotation] = field(default_factory=list, repr=False, compare=False)

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
            raise Exception('This VariantCall object does not encode a translocation. Therefore translocation orientation cannot be inferred.')

    def to_dict(self) -> Dict:
        data = {
            'variant_call_id': ['' if self.id is None else self.id],
            'source_id': ['' if self.source_id is None else self.source_id],
            'sample_id': ['' if self.sample_id is None else self.sample_id],
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
            'variant_sequences': [';'.join([i.sequence for i in self.variant_sequences])],
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

        pos_1_annotation_chrom = []
        pos_1_annotation_pos = []
        pos_1_annotation_region = []
        pos_1_annotation_gene_id = []
        pos_1_annotation_gene_source = []
        pos_1_annotation_gene_name = []
        pos_1_annotation_gene_chromosome = []
        pos_1_annotation_gene_start = []
        pos_1_annotation_gene_end = []
        pos_1_annotation_gene_strand = []
        pos_1_annotation_gene_type = []
        pos_1_annotation_gene_level = []
        pos_1_annotation_gene_version = []
        for i in self.position_1_annotations:
            pos_1_annotation_chrom.append(str(i.chrom))
            pos_1_annotation_pos.append(str(i.pos))
            pos_1_annotation_region.append(str(i.region))
            pos_1_annotation_gene_id.append('' if i.gene is None else str(i.gene.id))
            pos_1_annotation_gene_source.append('' if i.gene is None else str(i.gene.source))
            pos_1_annotation_gene_name.append('' if i.gene is None else str(i.gene.name))
            pos_1_annotation_gene_chromosome.append('' if i.gene is None else str(i.gene.chromosome))
            pos_1_annotation_gene_start.append('' if i.gene is None else str(i.gene.start))
            pos_1_annotation_gene_end.append('' if i.gene is None else str(i.gene.end))
            pos_1_annotation_gene_strand.append('' if i.gene is None else str(i.gene.strand))
            pos_1_annotation_gene_type.append('' if i.gene is None else str(i.gene.type))
            pos_1_annotation_gene_level.append('' if i.gene is None else str(i.gene.level))
            pos_1_annotation_gene_version.append('' if i.gene is None else str(i.gene.version))

        pos_2_annotation_chrom = []
        pos_2_annotation_pos = []
        pos_2_annotation_region = []
        pos_2_annotation_gene_id = []
        pos_2_annotation_gene_source = []
        pos_2_annotation_gene_name = []
        pos_2_annotation_gene_chromosome = []
        pos_2_annotation_gene_start = []
        pos_2_annotation_gene_end = []
        pos_2_annotation_gene_strand = []
        pos_2_annotation_gene_type = []
        pos_2_annotation_gene_level = []
        pos_2_annotation_gene_version = []
        for i in self.position_2_annotations:
            pos_2_annotation_chrom.append(str(i.chrom))
            pos_2_annotation_pos.append(str(i.pos))
            pos_2_annotation_region.append(str(i.region))
            pos_2_annotation_gene_id.append('' if i.gene is None else str(i.gene.id))
            pos_2_annotation_gene_source.append('' if i.gene is None else str(i.gene.source))
            pos_2_annotation_gene_name.append('' if i.gene is None else str(i.gene.name))
            pos_2_annotation_gene_chromosome.append('' if i.gene is None else str(i.gene.chromosome))
            pos_2_annotation_gene_start.append('' if i.gene is None else str(i.gene.start))
            pos_2_annotation_gene_end.append('' if i.gene is None else str(i.gene.end))
            pos_2_annotation_gene_strand.append('' if i.gene is None else str(i.gene.strand))
            pos_2_annotation_gene_type.append('' if i.gene is None else str(i.gene.type))
            pos_2_annotation_gene_level.append('' if i.gene is None else str(i.gene.level))
            pos_2_annotation_gene_version.append('' if i.gene is None else str(i.gene.version))

        data['position_1_annotation_chromosome'] = [';'.join(pos_1_annotation_chrom)]
        data['position_1_annotation_position'] = [';'.join(pos_1_annotation_pos)]
        data['position_1_annotation_region'] = [';'.join(pos_1_annotation_region)]
        data['position_1_annotation_gene_id'] = [';'.join(pos_1_annotation_gene_id)]
        data['position_1_annotation_gene_source'] = [';'.join(pos_1_annotation_gene_source)]
        data['position_1_annotation_gene_name'] = [';'.join(pos_1_annotation_gene_name)]
        data['position_1_annotation_gene_chromosome'] = [';'.join(pos_1_annotation_gene_chromosome)]
        data['position_1_annotation_gene_start'] = [';'.join(pos_1_annotation_gene_start)]
        data['position_1_annotation_gene_end'] = [';'.join(pos_1_annotation_gene_end)]
        data['position_1_annotation_gene_strand'] = [';'.join(pos_1_annotation_gene_strand)]
        data['position_1_annotation_gene_type'] = [';'.join(pos_1_annotation_gene_type)]
        data['position_1_annotation_gene_level'] = [';'.join(pos_1_annotation_gene_level)]
        data['position_1_annotation_gene_version'] = [';'.join(pos_1_annotation_gene_version)]

        data['position_2_annotation_chrom'] = [';'.join(pos_2_annotation_chrom)]
        data['position_2_annotation_pos'] = [';'.join(pos_2_annotation_pos)]
        data['position_2_annotation_region'] = [';'.join(pos_2_annotation_region)]
        data['position_2_annotation_gene_id'] = [';'.join(pos_2_annotation_gene_id)]
        data['position_2_annotation_gene_source'] = [';'.join(pos_2_annotation_gene_source)]
        data['position_2_annotation_gene_name'] = [';'.join(pos_2_annotation_gene_name)]
        data['position_2_annotation_gene_chromosome'] = [';'.join(pos_2_annotation_gene_chromosome)]
        data['position_2_annotation_gene_start'] = [';'.join(pos_2_annotation_gene_start)]
        data['position_2_annotation_gene_end'] = [';'.join(pos_2_annotation_gene_end)]
        data['position_2_annotation_gene_strand'] = [';'.join(pos_2_annotation_gene_strand)]
        data['position_2_annotation_gene_type'] = [';'.join(pos_2_annotation_gene_type)]
        data['position_2_annotation_gene_level'] = [';'.join(pos_2_annotation_gene_level)]
        data['position_2_annotation_gene_version'] = [';'.join(pos_2_annotation_gene_version)]
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())
