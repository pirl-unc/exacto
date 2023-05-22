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
The purpose of this python3 script is to implement the Gene dataclass.
"""


import pandas as pd
from bisect import bisect_left, bisect_right, insort
from dataclasses import dataclass, field
from typing import List, Dict
from .exon import Exon
from .logging import get_logger
from .transcript import Transcript


logger = get_logger(__name__)


@dataclass(frozen=True)
class Gene:
    id: str                                                         # e.g. 'ENSG00000122674.10'
    stable_id: str                                                  # e.g. 'ENSG00000122674'
    source: str                                                     # e.g. 'ENSEMBL'
    source_version: str                                             # e.g. '107'
    name: str                                                       # e.g. 'CCZ1'
    chromosome: str                                                 # e.g. 'chr7'
    start: int                                                      # e.g. 1000
    end: int                                                        # e.g. 10000
    strand: str                                                     # e.g. '+' or '-' (one of the Strands values in constants.py)
    type: str                                                       # e.g. 'protein_coding'
    level: str = None                                               # e.g. 1
    version: int = None                                             # e.g. 10
    genome: str = None                                              # e.g. 'GRCh38'
    transcripts: List[Transcript] = field(default_factory=list)

    def __lt__(self, other):
        if isinstance(other, Gene):
            return (self.chromosome, self.start) < (other.chromosome, other.start)
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Gene):
            return self.id == self.id
        return NotImplemented

    @property
    def transcripts_count(self):
        return len(self.transcripts)

    @property
    def transcript_ids(self) -> List[str]:
        return [transcript.id for transcript in self.transcripts]

    def add_transcript(self, transcript: Transcript):
        """
        Adds a Transcript object.

        Parameters
        ----------
        transcript  :   Transcript object.
        """
        insort(self.transcripts, transcript)

    def to_dict(self):
        data = {
            'gene_id': [self.id],
            'gene_stable_id': [self.stable_id],
            'gene_source': [self.source],
            'gene_source_version': [self.source_version],
            'gene_name': [self.name],
            'gene_chromosome': [self.chromosome],
            'gene_start': [self.start],
            'gene_end': [self.end],
            'gene_strand': [self.strand],
            'gene_type': [self.type],
            'gene_level': [self.level],
            'gene_version': [self.version],
            'gene_genome': [self.genome],
            'gene_transcripts_count': [self.transcripts_count]
        }
        return data

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.to_dict())
