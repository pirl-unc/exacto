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
The purpose of this python3 script is to implement the GeneSet class.
"""


import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, ClassVar
from .exon import Exon
from .transcript import Transcript
from .gene import Gene
from .logging import get_logger


logger = get_logger(__name__)


@dataclass
class GeneSet:
    genes: Dict = field(default_factory=dict) # key = gene ID, value = an instance of the Gene class

    @property
    def gene_ids(self) -> List[str]:
        return self.genes.keys()

    def add_gene(
            self,
            gene: Gene):
        """
        Adds a gene.

        Parameters
        ----------
        gene    :   An instance of the Gene class.
        """
        if gene.id not in self.genes.keys():
            self.genes[gene.id] = gene
        else:
            logger.error('Gene with ID %s already exists.' % gene.id)
            exit(1)

    def add_transcript(
            self,
            gene_id: str,
            transcript: Transcript):
        """
        Add a transcript.

        Parameters
        ----------
        gene_id     :   Gene ID.
        transcript  :   An instance of the Transcript class.
        """
        if gene_id in self.genes.keys():
            self.genes[gene_id].add_transcript(transcript=transcript)
        else:
            logger.error('Gene with ID %s does not exist.')
            exit(1)

    def add_exon(
            self,
            gene_id: str,
            transcript_id: str,
            exon: Exon):
        """
        Add an exon.

        Parameters
        ----------
        gene_id         :   Gene ID.
        transcript_id   :   Transcript ID.
        exon            :   An instance of the Exon class.
        """
        if gene_id in self.genes.keys():
            if transcript_id in self.genes[gene_id].transcript_ids:
                self.genes[gene_id].add_exon(exon=exon)
            else:
                logger.error('Transcript with ID %s does not exist.' % transcript_id)
        else:
            logger.error('Gene with ID %s does not exist.' % gene_id)
            exit(1)
