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
The purpose of this python3 script is to implement Variant dataclass.
"""


import statistics
import pandas as pd
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from bisect import bisect_left, bisect_right, insort
from .variant_annotation import VariantAnnotation
from .variant_call import VariantCall


@dataclass(frozen=True)
class Variant:
    id: str
    variant_calls: List[VariantCall] = field(default_factory=list)

    @property
    def num_variant_calls(self) -> int:
        return len(self.variant_calls)

    @property
    def variant_call_ids(self) -> List[str]:
        return [variant_call.id for variant_call in self.variant_calls]

    @property
    def source_ids(self) -> List[str]:
        return [i.source_id for i in self.variant_calls]

    @property
    def sample_ids(self) -> str:
        return [i.sample_id for i in self.variant_calls]

    @property
    def nucleic_acids(self) -> List[str]:
        return [i.nucleic_acid for i in self.variant_calls]

    @property
    def variant_calling_methods(self) -> List[str]:
        return [i.variant_calling_method for i in self.variant_calls]

    @property
    def sequencing_platform(self) -> List[str]:
        return [i.sequencing_platform for i in self.variant_calls]

    @property
    def chromosome_1(self):
        return self.variant_calls[0].chromosome_1

    @property
    def chromosome_2(self):
        return self.variant_calls[0].chromosome_2

    @property
    def position_1(self) -> List[int]:
        return [i.position_1 for i in self.variant_calls]

    @property
    def position_1_stdev(self) -> float:
        return 0.0 if len(self.position_1) == 1 else statistics.stdev(self.position_1)

    @property
    def position_2(self) -> List[int]:
        return [i.position_2 for i in self.variant_calls]

    @property
    def position_2_stdev(self) -> float:
        return 0.0 if len(self.position_2) == 1 else statistics.stdev(self.position_2)

    @property
    def reference_alleles(self) -> List[str]:
        return [i.reference_allele for i in self.variant_calls]

    @property
    def alternate_alleles(self) -> List[str]:
        return [i.alternate_allele for i in self.variant_calls]

    @property
    def filter(self) -> List[str]:
        return [i.filter for i in self.variant_calls]

    @property
    def quality_scores(self) -> List[float]:
        return [i.quality_score for i in self.variant_calls]

    @property
    def precise(self) -> List[bool]:
        return [i.precise for i in self.variant_calls]

    @property
    def variant_type(self) -> str:
        return self.variant_calls[0].variant_type

    @property
    def variant_subtypes(self) -> List[str]:
        return [i.variant_subtype for i in self.variant_calls]

    @property
    def variant_sizes(self) -> List[int]:
        return [i.variant_size for i in self.variant_calls]

    @property
    def variant_size_stdev(self) -> float:
        return 0.0 if len(self.variant_sizes) == 1 else statistics.stdev(self.variant_sizes)

    @property
    def variant_sequences(self) -> List[List[str]]:
        return [i.variant_sequence for i in self.variant_calls]

    @property
    def total_read_counts(self) -> List[int]:
        return [i.total_read_count for i in self.variant_calls]

    @property
    def reference_allele_read_counts(self) -> List[int]:
        return [i.reference_allele_read_count for i in self.variant_calls]

    @property
    def alternate_allele_read_counts(self) -> List[int]:
        return [i.alternate_allele_read_count for i in self.variant_calls]

    @property
    def alternate_allele_fractions(self) -> List[float]:
        return [i.alternate_allele_fraction for i in self.variant_calls]

    @property
    def alternate_allele_read_ids(self) -> List[List[str]]:
        return [i.alternate_allele_read_ids for i in self.variant_calls]

    @property
    def tool_attributes(self) -> List[OrderedDict]:
        return [i.tool_attributes for i in self.variant_calls]

    @property
    def position_1_annotations(self) -> List[List[VariantAnnotation]]:
        return [i.position_1_annotations for i in self.variant_calls]

    @property
    def position_2_annotations(self) -> List[List[VariantAnnotation]]:
        return [i.position_2_annotations for i in self.variant_calls]

    def __eq__(self, other):
        if isinstance(other, Variant):
            return (self.id) == (other.id)
        return NotImplemented

    def add_variant_call(self, variant_call: VariantCall):
        """
        Adds a VariantCall object.

        Parameters
        ----------
        variant_call    :   VariantCall object.
        """
        if len(self.variant_calls) > 0:
            if variant_call.chromosome_1 != self.variant_calls[0].chromosome_1 or \
                variant_call.chromosome_2 != self.variant_calls[0].chromosome_2:
                raise Exception('VariantCall chromosome_1 and chromosome_2 must be the same as '
                                'existing VariantCall objects in this Variant object.')
        insort(self.variant_calls, variant_call)

    def to_dict(self) -> Dict:
        data = {
            'id': self.id,
            'variant_calls': [variant_call.to_dict() for variant_call in self.variant_calls]
        }
        return data

    def to_dataframe_row(self) -> Dict:
        data = defaultdict(list)
        for variant_call in self.variant_calls:
            data['variant_id'].append(self.id)
            for key, value in variant_call.to_dataframe_row().items():
                data[key].append(value[0])
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())
