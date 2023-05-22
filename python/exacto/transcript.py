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
The purpose of this python3 script is to implement the Transcript dataclass.
"""


import pandas as pd
from bisect import insort
from dataclasses import dataclass, field
from typing import List
from .exon import Exon
from .nucleotide_sequence import NucleotideSequence


@dataclass(frozen=True)
class Transcript:
    id: str                                                     # e.g. 'ENST00000445840.10'
    stable_id: str                                              # e.g. 'ENST00000445840'
    source: str                                                 # e.g. 'ENSEMBL'
    source_version: str                                         # e.g. '107'
    chromosome: str                                             # e.g. 'chr7'
    start: int                                                  # e.g. 1000
    end: int                                                    # e.g. 10000
    type: str                                                   # e.g. 'protein_coding'
    strand: str                                                 # e.g. '+' or '-' (one of the Strands values in constants.py)
    transcription_start_site: int = None                        # e.g. 1000
    version: int = None                                         # e.g. 10
    name: str = None                                            # e.g. 'CCZ1-201'
    level: int = None                                           # e.g. 1
    support_level: int = None                                   # e.g. 1
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
    genome: str = None                                          # e.g. 'GRCh38'
    tags: List[str] = field(default_factory=list)
    exons: List[Exon] = field(default_factory=list)

    def __lt__(self, other):
        if isinstance(other, Transcript):
            return (self.start, self.end) < (other.start, other.end)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Transcript):
            return self.id == self.id
        return NotImplemented

    @property
    def cds_length(self) -> int:
        """
        Returns transcript CDS length.
        """
        length = 0
        for curr_exon in self.exons:
            length += curr_exon.length
        return length

    @property
    def exons_count(self) -> int:
        return len(self.exons)

    @property
    def exon_ids(self) -> List[str]:
        return [exon.id for exon in self.exons]

    @property
    def length(self) -> int:
        """
        Returns transcript length (including UTRs and CDS).
        """
        cds_length = self.cds_length
        cds_length += abs(self.five_prime_utr_start - self.five_prime_utr_end) + 1
        cds_length += abs(self.three_prime_utr_start - self.three_prime_utr_end) + 1
        return cds_length

    @property
    def sequence(self) -> NucleotideSequence:
        sequence = ''
        for exon in self.exons:
            sequence += exon.sequence.sequence
        return NucleotideSequence(sequence=sequence)

    def add_exon(self, exon: Exon):
        """
        Adds an exon.

        Parameters
        ----------
        exon        :   An instance of the Exon class.
        """
        insort(self.exons, exon)

    def to_dict(self):
        data = {
            'transcript_id': [self.id],
            'transcript_stable_id': [self.stable_id],
            'transcript_source': [self.source],
            'transcript_source_version': [self.source_version],
            'transcript_chromosome': [self.chromosome],
            'transcript_start': [self.start],
            'transcript_end': [self.end],
            'transcript_type': [self.type],
            'transcript_strand': [self.strand],
            'transcript_transcription_start_site': [self.transcription_start_site],
            'transcript_version': [self.version],
            'transcript_name': [self.name],
            'transcript_level': [self.level],
            'transcript_support_level': [self.support_level],
            'transcript_five_prime_utr_start': [self.five_prime_utr_start],
            'transcript_five_prime_utr_end': [self.five_prime_utr_end],
            'transcript_three_prime_utr_start': [self.three_prime_utr_start],
            'transcript_three_prime_utr_end': [self.three_prime_utr_end],
            'transcript_start_codon_start': [self.start_codon_start],
            'transcript_start_codon_end': [self.start_codon_end],
            'transcript_stop_codon_start': [self.stop_codon_start],
            'transcript_stop_codon_end': [self.stop_codon_end],
            'transcript_utr_start': [self.utr_start],
            'transcript_utr_end': [self.utr_end],
            'transcript_genome': [self.genome],
            'transcript_exons_count': [self.exons_count],
            'transcript_tags': [';'.join(self.tags)]
        }
        return data

    def to_dataframe(self):
        return pd.DataFrame(self.to_dict())

