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


from dataclasses import dataclass, field
from typing import List
from .exon import Exon
from .edit import Edit


@dataclass
class Match(Exon):
    gene_id: str = ''
    transcript_id: str = ''
    exon_id: str = ''
    exon_number: int = -1
    chrom: str = ''
    start: int = -1
    end: int = -1
    strand: str = ''
    length: int = -1
    sequence: str = ''
    edits: List = field(default_factory=lambda: [])

    def __init__(self,
                 gene_id: str,
                 transcript_id: str,
                 exon_id: str,
                 exon_number: int,
                 strand: str,
                 chrom: str,
                 start: int,
                 end: int,
                 length: int,
                 sequence: str):
        super().__init__()
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.exon_id = exon_id
        self.exon_number = exon_number
        self.strand = strand
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = length
        self.sequence = sequence
        self.edits = []

        # Append edits
        positions = range(self.start, self.end + 1)
        for i in range(0, len(positions)):
            edit = Edit(
                ref=self.sequence[i],
                alt=self.sequence[i],
                pos=positions[i],
                sequence=self.sequence[i]
            )
            self.edits.append(edit)

    def __str__(self):
        msg = "[REFERENCE]\n"
        msg += "\tgene ID\t\t\t:\t%s\n" % self.gene_id
        msg += "\ttranscript ID\t\t:\t%s\n" % self.transcript_id
        msg += "\texon ID\t\t\t:\t%s\n" % self.exon_id
        msg += "\texon number\t\t:\t%i\n" % self.exon_number
        msg += "\tstrand\t\t\t:\t%s\n" % self.strand
        msg += "\tchromosome\t\t:\t%s\n" % self.chrom
        msg += "\tstart\t\t\t:\t%i\n" % self.start
        msg += "\tend\t\t\t:\t%i\n" % self.end
        msg += "\tlength\t\t\t:\t%i\n" % self.length
        msg += "\treference sequence\t:\t%s\n" % self.sequence
        return msg + \
               super(Match, self).__str__()
