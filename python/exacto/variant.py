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
The purpose of this python3 script is to implement Variant class.
"""


import statistics
import pandas as pd
from collections import defaultdict
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .variant_call import VariantCall
from .variant_annotation import VariantAnnotation


@dataclass
class Variant:
    id: str = None
    variant_calls: List[VariantCall] = field(default_factory=list)

    @property
    def variant_call_id(self) -> List[str]:
        return [i.id for i in self.variant_calls]

    @property
    def source_id(self) -> str:
        return [i.source_id for i in self.variant_calls]

    @property
    def tumor_sample_id(self) -> str:
        return [i.tumor_sample_id for i in self.variant_calls]

    @property
    def normal_sample_id(self) -> List[str]:
        return [i.normal_sample_id for i in self.variant_calls]

    @property
    def nucleic_acid(self) -> List[str]:
        return [i.nucleic_acid for i in self.nucleic_acid]

    @property
    def variant_calling_method(self) -> List[str]:
        return [i.variant_calling_method for i in self.variant_calls]

    @property
    def sequencing_platform(self) -> List[str]:
        return [i.sequencing_platform for i in self.variant_calls]

    @property
    def chr_1(self) -> str:
        return self.variant_calls[0].chr_1

    @property
    def pos_1(self) -> List[int]:
        return [i.pos_1 for i in self.variant_calls]

    @property
    def pos_1_stdev(self) -> float:
        return 0.0 if len(self.pos_1) == 1 else statistics.stdev(self.pos_1)

    @property
    def chr_2(self) -> str:
        return self.variant_calls[0].chr_2

    @property
    def pos_2(self) -> List[int]:
        return [i.pos_2 for i in self.variant_calls]

    @property
    def pos_2_stdev(self) -> float:
        return 0.0 if len(self.pos_2) == 1 else statistics.stdev(self.pos_2)

    @property
    def ref(self) -> List[str]:
        return [i.ref for i in self.variant_calls]

    @property
    def alt(self) -> List[str]:
        return [i.alt for i in self.variant_calls]

    @property
    def filter(self) -> List[str]:
        return [i.filter for i in self.variant_calls]

    @property
    def quality_score(self) -> List[float]:
        return [i.quality_score for i in self.variant_calls]

    @property
    def precise(self) -> List[bool]:
        return [i.precise for i in self.variant_calls]

    @property
    def variant_type(self) -> List[str]:
        return [i.variant_type for i in self.variant_calls]

    @property
    def variant_subtype(self) -> List[str]:
        return [i.variant_subtype for i in self.variant_calls]

    @property
    def variant_size(self) -> List[int]:
        return [i.variant_size for i in self.variant_calls]

    @property
    def variant_sequence(self) -> List[List[str]]:
        return [i.variant_sequence for i in self.variant_calls]

    @property
    def total_tumor_reads(self) -> List[int]:
        return [i.total_tumor_reads for i in self.variant_calls]

    @property
    def ref_tumor_reads(self) -> List[int]:
        return [i.ref_tumor_reads for i in self.variant_calls]

    @property
    def alt_tumor_reads(self) -> List[int]:
        return [i.alt_tumor_reads for i in self.variant_calls]

    @property
    def other_tumor_reads(self) -> List[int]:
        return [i.other_tumor_reads for i in self.variant_calls]

    @property
    def alt_tumor_fraction(self) -> List[float]:
        return [i.alt_tumor_fraction for i in self.variant_calls]

    @property
    def total_normal_reads(self) -> List[int]:
        return [i.total_normal_reads for i in self.variant_calls]

    @property
    def ref_normal_reads(self) -> List[int]:
        return [i.ref_normal_reads for i in self.variant_calls]

    @property
    def alt_normal_reads(self) -> List[int]:
        return [i.alt_normal_reads for i in self.variant_calls]

    @property
    def other_normal_reads(self) -> List[int]:
        return [i.other_normal_reads for i in self.variant_calls]

    @property
    def alt_normal_fraction(self) -> List[float]:
        return [i.alt_normal_fraction for i in self.variant_calls]

    @property
    def alt_tumor_read_id(self) -> List[List[str]]:
        return [i.alt_tumor_read_id for i in self.variant_calls]

    @property
    def alt_normal_read_id(self) -> List[List[str]]:
        return [i.alt_normal_read_id for i in self.variant_calls]

    @property
    def alt_tumor_softclip_direction(self) -> List[str]:
        return [i.alt_tumor_softclip_direction for i in self.variant_calls]

    @property
    def alt_normal_softclip_direction(self) -> List[str]:
        return [i.alt_normal_softclip_direction for i in self.variant_calls]

    @property
    def tool_attributes(self) -> List[OrderedDict]:
        return [i.tool_attributes for i in self.variant_calls]

    @property
    def pos_1_annotations(self) -> List[List[VariantAnnotation]]:
        return [i.pos_1_annotations for i in self.variant_calls]

    @property
    def pos_2_annotations(self) -> List[List[VariantAnnotation]]:
        return [i.pos_2_annotations for i in self.variant_calls]

    def to_dict(self) -> Dict:
        data = defaultdict(list)
        for variant_call in self.variant_calls:
            data['variant_id'].append(self.id)
            for key, value in variant_call.to_dict().items():
                data[key].append(value[0])
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())
