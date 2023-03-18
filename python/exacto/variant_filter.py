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
The purpose of this python3 script is to implement the Filter dataclass.
"""


import numpy as np
from dataclasses import dataclass, field
from .constants import VariantFilterQuantifiers, VariantFilterOperators
from .variant import Variant
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class VariantFilter:
    quantifier: str         # 'all', 'any', 'min', 'max', median', 'average'
    attribute: str          # 'alt_tumor_reads'
    operator: str           # '<', '<=', '>', '>=', '==', 'in'
    value: None             # 3, ["chr1","chr2","chr3"]

    def is_predicate(self, variant: Variant) -> bool:
        """
        Applies filter to the input variant and returns True if the variant
        meets the filter criterion and returns False if the variant does not
        meet the filter criterion.

        Parameters
        ----------
        variant     :   An instance of the Variant class.

        Returns
        -------
        True or False.
        """
        if self.quantifier == VariantFilterQuantifiers.ALL:
            query = '%s %s %s' % (self.attribute, self.operator, self.value)
            for variant_call in variant.variant_calls:
                df_variant_call = variant_call.to_dataframe()
                df_variant_call_match = df_variant_call.query(query)
                if len(df_variant_call_match) == 0:
                    return False
            return True
        elif self.quantifier == VariantFilterQuantifiers.ANY:
            query = '%s %s %s' % (self.attribute, self.operator, self.value)
            for variant_call in variant.variant_calls:
                df_variant_call = variant_call.to_dataframe()
                df_variant_call_match = df_variant_call.query(query)
                if len(df_variant_call_match) > 0:
                    return True
            return False
        elif self.quantifier == VariantFilterQuantifiers.AVERAGE:
            attribute_values = []
            for variant_call in variant.variant_calls:
                df_variant_call = variant_call.to_dataframe()
                attribute_values.append(df_variant_call[self.attribute].values.tolist()[0])
            summarized_value = np.mean(attribute_values)
            if self.operator == VariantFilterOperators.LESS_THAN:
                return True if summarized_value < self.value else False
            elif self.operator == VariantFilterOperators.LESS_THAN_OR_EQUAL_TO:
                return True if summarized_value <= self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN:
                return True if summarized_value > self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO:
                return True if summarized_value >= self.value else False
            else:
                logger.error('Unknown operator: %s.' % self.operator)
                exit(1)
        elif self.quantifier == VariantFilterQuantifiers.MEDIAN:
            attribute_values = []
            for variant_call in variant.variant_calls:
                df_variant_call = variant_call.to_dataframe()
                attribute_values.append(df_variant_call[self.attribute].values.tolist()[0])
            summarized_value = np.median(attribute_values)
            if self.operator == VariantFilterOperators.LESS_THAN:
                return True if summarized_value < self.value else False
            elif self.operator == VariantFilterOperators.LESS_THAN_OR_EQUAL_TO:
                return True if summarized_value <= self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN:
                return True if summarized_value > self.value else False
            elif self.operator == VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO:
                return True if summarized_value >= self.value else False
            else:
                logger.error('Unknown operator: %s.' % self.operator)
                exit(1)
        else:
            logger.error('Unknown quantifier: %s.' % self.quantifier)
            exit(1)
