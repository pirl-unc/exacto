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
The purpose of this python3 script is to implement the GenomicRange dataclass.
"""


import pandas as pd
from dataclasses import dataclass, field
from functools import total_ordering
from typing import Dict
from .logging import get_logger


logger = get_logger(__name__)


@total_ordering
@dataclass(frozen=True)
class GenomicRange:
    chromosome: str
    start: int
    end: int

    @property
    def id(self):
        return '%s:%i-%i' % (self.chromosome, self.start, self.end)

    def __lt__(self, other):
        if isinstance(other, GenomicRange):
            return (self.chromosome, self.start, self.end) < \
                   (other.chromosome, other.start, self.end)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, GenomicRange):
            return (self.chromosome, self.start, self.end) == \
                   (other.chromosome, other.start, self.end)
        return NotImplemented

    def overlaps(self, chromosome: str, start: int, end: int) -> bool:
        """
        Returns True if query position overlaps the GenomicRange.

        Parameters
        ----------
        chromosome  :   Chromosome.
        start       :   Start position.
        end         :   End position.

        Returns
        -------
        True or False.
        """
        # By De Morgan's law on checking for non-overlapping regions
        if chromosome == self.chromosome and \
            start <= self.end and \
            end >= self.start:
            return True
        else:
            return False

    def to_dict(self) -> Dict:
        data = {
            'chromosome': self.chromosome,
            'start': self.start,
            'end': self.end
        }
        return data

    def to_dataframe_row(self) -> Dict:
        data = {
            'chromosome': [self.chromosome],
            'start': [self.start],
            'end': [self.end]
        }
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dataframe_row())
