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
The purpose of this python3 script is to implement the VariantFilter dataclass.
"""


import json
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Type, Dict
from .constants import VariantFilterQuantifiers, VariantFilterOperators
from .logging import get_logger
from .variant import Variant


logger = get_logger(__name__)


@dataclass
class VariantFilter:
    quantifier: str                                         # 'all', 'any', 'min', 'max', median', 'average'
    attribute: str                                          # 'alternate_allele_read_count'
    operator: str                                           # '<', '<=', '>', '>=', '==', '!=', 'in'
    value: None                                             # 3, ["chr1","chr2","chr3"]
    sample_ids: List[str] = field(default_factory=list)     # sample IDs

    def __post_init__(self):
        # Check if value is a list
        if type(self.value) == str:
            if self.value[0] == '[' and self.value[-1] == ']':
                parsed_list = json.loads(self.value)
                if isinstance(parsed_list, list):
                    self.value = parsed_list

    def to_dict(self):
        data = {
            'quantifier': self.quantifier,
            'attribute': self.attribute,
            'operator': self.operator,
            'value': self.value,
            'sample_ids': self.sample_ids
        }
        return data

    def to_dataframe_row(self) -> Dict:
        data = {
            'quantifier': [self.quantifier],
            'attribute': [self.attribute],
            'operator': [self.operator],
            'value': [str(self.value)],
            'sample_ids': [';'.join(self.sample_ids)]
        }
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())
