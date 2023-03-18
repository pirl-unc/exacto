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


import pandas as pd
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .variant_annotation import VariantAnnotation


@dataclass
class VariantCall:
    id: str = None
    source_id: str = None
    tumor_sample_id: str = None
    normal_sample_id: str = None
    nucleic_acid: str = None
    variant_calling_method: str = None
    sequencing_platform: str = None
    chr_1: str = None
    pos_1: int = None
    chr_2: str = None
    pos_2: int = None
    ref: str = None
    alt: str = None
    filter: str = None
    quality_score: float = None
    precise: bool = None
    variant_type: str = None
    variant_subtype: str = None
    variant_size: int = None
    variant_sequences: List[str] = field(default_factory=list)
    total_tumor_reads: int = None
    ref_tumor_reads: int = None
    alt_tumor_reads: int = None
    other_tumor_reads: int = None
    alt_tumor_fraction: float = None
    total_normal_reads: int = None
    ref_normal_reads: int = None
    alt_normal_reads: int = None
    other_normal_reads: int = None
    alt_normal_fraction: float = None
    alt_tumor_read_ids: List[str] = field(default_factory=list)
    alt_normal_read_ids: List[str] = field(default_factory=list)
    alt_tumor_softclip_direction: str = None                                    # ++, +-, -+, --    <pos_1,pos_2>
    alt_normal_softclip_direction: str = None                                   # ++, +-, -+, --    <pos_1,pos_2>
    tool_attributes: OrderedDict = field(default_factory=dict)
    pos_1_annotations: List[VariantAnnotation] = field(default_factory=list)
    pos_2_annotations: List[VariantAnnotation] = field(default_factory=list)

    def to_dict(self) -> Dict:
        data = {
            'variant_call_id': ['' if self.id is None else self.id],
            'source_id': ['' if self.source_id is None else self.source_id],
            'tumor_sample_id': ['' if self.tumor_sample_id is None else self.tumor_sample_id],
            'normal_sample_id': ['' if self.normal_sample_id is None else self.normal_sample_id],
            'nucleic_acid': ['' if self.nucleic_acid is None else self.nucleic_acid],
            'variant_calling_method': ['' if self.variant_calling_method is None else self.variant_calling_method],
            'sequencing_platform': ['' if self.sequencing_platform is None else self.sequencing_platform],
            'chr_1': ['' if self.chr_1 is None else self.chr_1],
            'pos_1': ['' if self.pos_1 is None else self.pos_1],
            'chr_2': ['' if self.chr_2 is None else self.chr_2],
            'pos_2': ['' if self.pos_2 is None else self.pos_2],
            'ref': ['' if self.ref is None else self.ref],
            'alt': ['' if self.alt is None else self.alt],
            'filter': ['' if self.filter is None else self.filter],
            'quality_score': ['' if self.quality_score is None else self.quality_score],
            'precise': ['' if self.precise is None else self.precise],
            'variant_type': ['' if self.variant_type is None else self.variant_type],
            'variant_subtype': ['' if self.variant_subtype is None else self.variant_subtype],
            'variant_size': ['' if self.variant_size is None else self.variant_size],
            'variant_sequences': [';'.join(self.variant_sequences)],
            'total_tumor_reads': ['' if self.total_tumor_reads is None else self.total_tumor_reads],
            'ref_tumor_reads': ['' if self.ref_tumor_reads is None else self.ref_tumor_reads],
            'alt_tumor_reads': ['' if self.alt_tumor_reads is None else self.alt_tumor_reads],
            'other_tumor_reads': ['' if self.other_tumor_reads is None else self.other_tumor_reads],
            'alt_tumor_fraction': ['' if self.alt_tumor_fraction is None else self.alt_tumor_fraction],
            'total_normal_reads': ['' if self.total_normal_reads is None else self.total_normal_reads],
            'ref_normal_reads': ['' if self.ref_normal_reads is None else self.ref_normal_reads],
            'alt_normal_reads': ['' if self.alt_normal_reads is None else self.alt_normal_reads],
            'other_normal_reads': ['' if self.other_normal_reads is None else self.other_normal_reads],
            'alt_normal_fraction': ['' if self.alt_normal_fraction is None else self.alt_normal_fraction],
            'alt_tumor_read_ids': [';'.join(self.alt_tumor_read_ids)],
            'alt_normal_read_ids': [';'.join(self.alt_normal_read_ids)],
            'alt_tumor_softclip_direction': ['' if self.alt_tumor_softclip_direction is None else self.alt_tumor_softclip_direction],
            'alt_normal_softclip_direction': ['' if self.alt_normal_softclip_direction is None else self.alt_normal_softclip_direction],
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
        for i in self.pos_1_annotations:
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
        for i in self.pos_2_annotations:
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

        data['pos_1_annotation_chrom'] = [';'.join(pos_1_annotation_chrom)]
        data['pos_1_annotation_pos'] = [';'.join(pos_1_annotation_pos)]
        data['pos_1_annotation_region'] = [';'.join(pos_1_annotation_region)]
        data['pos_1_annotation_gene_id'] = [';'.join(pos_1_annotation_gene_id)]
        data['pos_1_annotation_gene_source'] = [';'.join(pos_1_annotation_gene_source)]
        data['pos_1_annotation_gene_name'] = [';'.join(pos_1_annotation_gene_name)]
        data['pos_1_annotation_gene_chromosome'] = [';'.join(pos_1_annotation_gene_chromosome)]
        data['pos_1_annotation_gene_start'] = [';'.join(pos_1_annotation_gene_start)]
        data['pos_1_annotation_gene_end'] = [';'.join(pos_1_annotation_gene_end)]
        data['pos_1_annotation_gene_strand'] = [';'.join(pos_1_annotation_gene_strand)]
        data['pos_1_annotation_gene_type'] = [';'.join(pos_1_annotation_gene_type)]
        data['pos_1_annotation_gene_level'] = [';'.join(pos_1_annotation_gene_level)]
        data['pos_1_annotation_gene_version'] = [';'.join(pos_1_annotation_gene_version)]

        data['pos_2_annotation_chrom'] = [';'.join(pos_2_annotation_chrom)]
        data['pos_2_annotation_pos'] = [';'.join(pos_2_annotation_pos)]
        data['pos_2_annotation_region'] = [';'.join(pos_2_annotation_region)]
        data['pos_2_annotation_gene_id'] = [';'.join(pos_2_annotation_gene_id)]
        data['pos_2_annotation_gene_source'] = [';'.join(pos_2_annotation_gene_source)]
        data['pos_2_annotation_gene_name'] = [';'.join(pos_2_annotation_gene_name)]
        data['pos_2_annotation_gene_chromosome'] = [';'.join(pos_2_annotation_gene_chromosome)]
        data['pos_2_annotation_gene_start'] = [';'.join(pos_2_annotation_gene_start)]
        data['pos_2_annotation_gene_end'] = [';'.join(pos_2_annotation_gene_end)]
        data['pos_2_annotation_gene_strand'] = [';'.join(pos_2_annotation_gene_strand)]
        data['pos_2_annotation_gene_type'] = [';'.join(pos_2_annotation_gene_type)]
        data['pos_2_annotation_gene_level'] = [';'.join(pos_2_annotation_gene_level)]
        data['pos_2_annotation_gene_version'] = [';'.join(pos_2_annotation_gene_version)]
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())
