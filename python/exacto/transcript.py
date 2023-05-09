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
The purpose of this python3 script is to implement the Transcript class.
"""


from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar, Set
from bisect import insort
from .exon import Exon


@dataclass
class Transcript:
    id: str                                                     # e.g. 'ENST00000445840'
    chromosome: str                                             # e.g. 'chr7'
    start: int                                                  # e.g. 1000
    end: int                                                    # e.g. 10000
    type: str                                                   # e.g. 'protein_coding' (one of the TranscriptTypes values in constants.py)
    strand: str                                                 # e.g. '+' or '-' (one of the Strands values in constants.py)
    transcript_start_site: int = None                           # e.g. 1000
    source: str = None                                          # e.g. 'ENSEMBL'
    version: int = None                                         # e.g. 1
    name: str = None                                            # e.g. 'CCZ1-201'
    level: int = None                                           # e.g. 1
    five_prime_utr_start: int = None                            # e.g. 1000 (5' position)
    five_prime_utr_end: int = None                              # e.g. 1500 (3' position)
    three_prime_utr_start: int = None                           # e.g. 10000 (5' position)
    three_prime_utr_end: int = None                             # e.g. 10500 (3' position)
    start_codon_start: int = None                               # e.g. 1000
    start_codon_end: int = None                                 # e.g. 1200
    stop_codon_start: int = None                                # e.g. 9500
    stop_codon_end: int = None                                  # e.g. 10000
    utr_start: int = None                                       # e.g. 10000
    utr_end: int = None                                         # e.g. 10500
    _sequence: str = None                                       # e.g. 'ATA...CGT'
    exons: List[Exon] = field(default_factory=list)

    @property
    def sequence(self):
        if self._sequence is None:
            sequence = ''
            for exon in self.exons:
                sequence += exon.sequence
        else:
            return self._sequence

    @sequence.setter
    def sequence(self, sequence: str) -> None:
        self._sequence = sequence

    @property
    def exon_ids(self) -> List[str]:
        return [exon.id for exon in self.exons]

    @property
    def length(self):
        """
        Returns transcript length (including UTRs and CDS).
        """
        cds_length = self.cds_length
        cds_length += abs(self.five_prime_utr_start - self.five_prime_utr_end) + 1
        cds_length += abs(self.three_prime_utr_start - self.three_prime_utr_end) + 1
        return cds_length

    @property
    def cds_length(self):
        """
        Returns transcript CDS length.
        """
        length = 0
        for curr_exon in self.exons:
            length += curr_exon.length
        return length

    def add_exon(self, exon: Exon):
        """
        Adds an exon.

        Parameters
        ----------
        exon        :   An instance of the Exon class.
        """
        insort(self.exons, exon)
