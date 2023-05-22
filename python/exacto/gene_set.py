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


from bisect import bisect_left, bisect_right, insort
from dataclasses import dataclass, field
from typing import List, Dict
from .exon import Exon
from .gene import Gene
from .logging import get_logger
from .transcript import Transcript


logger = get_logger(__name__)


@dataclass
class GeneSet:
    genes: List[Gene] = field(default_factory=list)

    @property
    def gene_ids(self) -> List[str]:
        return [gene.id for gene in self.genes]

    def add_gene(self, gene: Gene):
        """
        Add a Gene object.

        Parameters
        ----------
        gene    :   Gene object.
        """
        insort(self.genes, gene)

    def add_transcript(self, gene_id: str, transcript: Transcript):
        """
        Add a Transcript object.

        Parameters
        ----------
        gene_id     :   Gene ID.
        transcript  :   An instance of the Transcript class.
        """
        index = self.genes.index(Gene(
            id=gene_id,
            stable_id=gene_id,
            source='',
            source_version='',
            name='',
            chromosome='',
            start=0,
            end=0,
            strand='',
            type='',
            level='',
            version=0,
            genome=''
        ))
        self.genes[index].add_transcript(transcript=transcript)

    def add_exon(self, gene_id: str, transcript_id: str, exon: Exon):
        """
        Add a Exon object.

        Parameters
        ----------
        gene_id         :   Gene ID.
        transcript_id   :   Transcript ID.
        exon            :   Exon object.
        """
        gene_index = self.genes.index(Gene(
            id=gene_id,
            stable_id=gene_id,
            source='',
            source_version='',
            name='',
            chromosome='',
            start=0,
            end=0,
            strand='',
            type='',
            level='',
            version=0,
            genome=''
        ))
        transcript_index = self.genes[gene_index].transcripts.index(Transcript(
            id=transcript_id,
            stable_id=transcript_id,
            source='',
            source_version='',
            chromosome='',
            start=0,
            end=0,
            type='',
            strand=''
        ))
        self.genes[gene_index].transcripts[transcript_index].add_exon(exon=exon)

    def find_genes(self, chromosome, position) -> List[Gene]:
        """
        Finds Gene objects that match the query parameters.

        Parameters
        ----------
        chromosome      :   Chromosome.
        position        :   Position.

        Returns
        -------
        genes           :   List of Gene objects.
        """
        # Find the leftmost index where chromosome, start, and end positions match
        left_index = bisect_left(
            self.genes,
            Gene(
                id='',
                stable_id='',
                source='',
                source_version='',
                name='',
                chromosome=chromosome,
                start=position,
                end=position,
                strand='',
                type='',
                level='',
                version='',
                genome=''
            )
        )

        # Find the rightmost index where chromosome_1 and chromosome_2 match
        right_index = bisect_right(
            self.genes,
            Gene(
                id='',
                stable_id='',
                source='',
                source_version='',
                name='',
                chromosome=chromosome,
                start=position,
                end=position,
                strand='',
                type='',
                level='',
                version='',
                genome=''
            )
        )

        genes = []
        for gene in self.genes[left_index:right_index]:
            if chromosome == gene.chromosome and gene.start <= position <= gene.end:
                genes.append(gene)
        return genes
