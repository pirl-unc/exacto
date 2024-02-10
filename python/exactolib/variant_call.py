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


from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class VariantCall:
    id: str
    chromosome_1: str
    position_1: int
    chromosome_2: str
    position_2: int
    variant_type: str
    reference_allele: str
    alternate_allele: str
    variant_size: int
    alternate_allele_read_ids: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        data = {
            'id': self.id,
            'chromosome_1': self.chromosome_1,
            'position_1': self.position_1,
            'chromosome_2': self.chromosome_2,
            'position_2': self.position_2,
            'variant_type': self.variant_type,
            'reference_allele': self.reference_allele,
            'alterante_allele': self.alternate_allele,
            'variant_size': self.variant_size,
            'alternate_allele_read_count': len(self.alternate_allele_read_ids),
            'alternate_allele_read_ids': ';'.join(self.alternate_allele_read_ids)
        }
        return data
