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
The purpose of this python3 script is to implement the Gene class.
"""


import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .exon import Exon
from .transcript import Transcript
from exacto.logging import get_logger


logger = get_logger(__name__)


@dataclass
class Gene:
    id: str                                             # e.g. 'ENSG00000122674'
    source: str                                         # e.g. 'ENSEMBL'
    name: str                                           # e.g. 'CCZ1'
    chromosome: str                                     # e.g. 'chr7'
    start: int                                          # e.g. 1000
    end: int                                            # e.g. 10000
    strand: str                                         # e.g. '+' or '-' (one of the Strands values in constants.py)
    type: str                                           # e.g. 'protein_coding' (one of the TranscriptTypes values in constants.py)
    level: str = None                                   # e.g. 1
    version: int = None                                 # e.g. 1
    transcripts: Dict = field(default_factory=dict)     # key = transcript ID, value = an instance of the Transcript class

    @property
    def transcript_ids(self):
        return [transcript.id for transcript in self.transcripts.values()]

    def add_transcript(
            self,
            transcript: Transcript):
        """
        Adds a transcript.

        Parameters
        ----------
        transcript  :   An instance of the Transcript class.
        """
        if transcript.id not in self.transcripts.keys():
            self.transcripts[transcript.id] = transcript
        else:
            logger.error('Transcript with ID %s already exists.' % transcript.id)
            exit(1)

    def add_exon(
            self,
            transcript_id: str,
            exon: Exon):
        if transcript_id in self.transcripts.keys():
            self.transcripts[transcript_id].add_exon(exon=exon)
        else:
            logger.error('Transcript with ID %s does not exist.' % transcript_id)
            exit(1)

    def to_dataframe(self) -> pd.DataFrame:
        df = pd.DataFrame({
            'gene_id': [self.id],
            'gene_source': [self.source],
            'gene_name': [self.name],
            'gene_chromosome': [self.chromosome],
            'gene_start': [self.start],
            'gene_end': [self.end],
            'gene_strand': [self.strand],
            'gene_type': [self.type],
            'gene_level': [self.level],
            'gene_version': [self.version]
        })
        return df
