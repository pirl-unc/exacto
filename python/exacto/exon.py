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
The purpose of this python3 script is to implement the Exon class.
"""


from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from functools import total_ordering


@total_ordering
@dataclass(frozen=True)
class Exon:
    id: str                                         # e.g. 'ENSE00001910415'
    chromosome: str                                 # e.g. 'chr1'
    start: int                                      # e.g. 1000 (5' start)
    end: int                                        # e.g. 1200 (3' end)
    sequence: str                                   # e.g. 'AAATCCT...TTGG' (5' to 3' sequence)
    number: int                                     # e.g. 1
    strand: str                                     # e.g. '+'
    version: int                                    # e.g. 1
    source: str = None                              # e.g. 'ENSEMBL'
    tags: List[str] = field(default_factory=list)   # e.g. 'Ensembl_canonical'

    @property
    def length(self) -> int:
        """
        Returns exon.
        """
        return len(self.sequence)

    def __lt__(self, other):
        if isinstance(other, Exon):
            return self.number < other.number
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Exon):
            return self.number == other.number
        return NotImplemented
