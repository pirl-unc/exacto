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
The purpose of this python3 script is to implement the Exon dataclass.
"""


import pandas as pd
from dataclasses import dataclass, field
from functools import total_ordering
from typing import List
from .nucleotide_sequence import NucleotideSequence


@total_ordering
@dataclass(frozen=True)
class Exon:
    id: str                                         # e.g. 'ENSE00001910415.10'
    stable_id: str                                  # e.g. 'ENSE00001910415'
    source: str                                     # e.g. 'ENSEMBL'
    source_version: str                             # e.g. '107'
    chromosome: str                                 # e.g. 'chr1'
    start: int                                      # e.g. 1000 (5' start)
    end: int                                        # e.g. 1200 (3' end)
    number: int                                     # e.g. 1
    strand: str                                     # e.g. '+'
    sequence: NucleotideSequence
    utr_start: int = None                           # e.g. 10000
    utr_end: int = None                             # e.g. 10500
    version: int = None                             # e.g. 10
    tags: List[str] = field(default_factory=list)   # e.g. 'Ensembl_canonical'

    @property
    def length(self) -> int:
        return len(self.sequence.sequence)

    def __lt__(self, other):
        if isinstance(other, Exon):
            return self.number < other.number
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Exon):
            return self.id == other.id
        return NotImplemented

    def to_dict(self):
        data = {
            'exon_id': [self.id],
            'exon_stable_id': [self.stable_id],
            'exon_source': [self.source],
            'exon_source_version': [self.source_version],
            'exon_chromosome': [self.chromosome],
            'exon_start': [self.start],
            'exon_end': [self.end],
            'exon_number': [self.number],
            'exon_strand': [self.strand],
            'exon_sequence': [self.sequence.sequence],
            'exon_utr_start': ['' if self.utr_start is None else self.utr_start],
            'exon_utr_end': ['' if self.utr_end is None else self.utr_end],
            'exon_version': [self.version],
            'exon_tags': [';'.join(self.tags)]
        }
        return data

    def to_dataframe(self):
        return pd.DataFrame(self.to_dict())

