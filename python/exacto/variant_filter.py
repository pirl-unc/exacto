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


import numpy as np
from dataclasses import dataclass, field
from typing import List
from .constants import VariantFilterQuantifiers, VariantFilterOperators
from .logging import get_logger
from .variant import Variant


logger = get_logger(__name__)


@dataclass
class VariantFilter:
    quantifier: str                                         # 'all', 'any', 'min', 'max', median', 'average'
    attribute: str                                          # 'alternate_allele_read_count'
    operator: str                                           # '<', '<=', '>', '>=', '==', 'in'
    value: None                                             # 3, ["chr1","chr2","chr3"]
    sample_ids: List[str] = field(default_factory=list)     # sample IDs

    def keep(self, variant: Variant) -> bool:
        """
        Applies filter to the input variant and returns True if the variant
        meets the filter criterion and returns False if the variant does not
        meet the filter criterion.

        Parameters
        ----------
        variant         :   Variant object.

        Returns
        -------
        True or False.
        """
        if self.quantifier == VariantFilterQuantifiers.ALL:
            query = '%s %s %s' % (self.attribute, self.operator, self.value)
            for variant_call in variant.variant_calls:
                if variant_call.sample_id in self.sample_ids:
                    try:
                        if len(variant_call.to_dataframe().query(query)) == 0:
                            return False
                    except:
                        return False
            return True
        elif self.quantifier == VariantFilterQuantifiers.ANY:
            query = '%s %s %s' % (self.attribute, self.operator, self.value)
            for variant_call in variant.variant_calls:
                if variant_call.sample_id in self.sample_ids:
                    try:
                        if len(variant_call.to_dataframe().query(query)) > 0:
                            return True
                    except:
                        pass
            return False
        else:
            attribute_values = []
            for variant_call in variant.variant_calls:
                if variant_call.sample_id in self.sample_ids:
                    df_variant_call = variant_call.to_dataframe()
                    attribute_values.append(df_variant_call[self.attribute].values.tolist()[0])

            if self.quantifier == VariantFilterQuantifiers.MIN:
                summarized_value = min(attribute_values)
            elif self.quantifier == VariantFilterQuantifiers.MAX:
                summarized_value = max(attribute_values)
            elif self.quantifier == VariantFilterQuantifiers.AVERAGE:
                summarized_value = np.mean(attribute_values)
            elif self.quantifier == VariantFilterQuantifiers.MEDIAN:
                summarized_value = np.median(attribute_values)
            else:
                raise Exception('Unknown quantifier: %s' % self.quantifier)

            if self.operator == VariantFilterOperators.LESS_THAN:
                return True if summarized_value < self.value else False
            elif self.operator == VariantFilterOperators.LESS_THAN_OR_EQUAL_TO:
                return True if summarized_value <= self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN:
                return True if summarized_value > self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO:
                return True if summarized_value >= self.value else False
            elif self.operator == VariantFilterOperators.EQUALS:
                return True if summarized_value == self.value else False
            else:
                raise Exception('Unknown operator: %s' % self.operator)
