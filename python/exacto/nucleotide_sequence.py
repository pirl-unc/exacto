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
The purpose of this python3 script is to implement the NucleotideSequence dataclass.
"""


from dataclasses import dataclass
from .constants import RNA_CODONS, AMINO_ACID_CODES
from .peptide_sequence import PeptideSequence


@dataclass
class NucleotideSequence:
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        gc_count = self.sequence.count('G') + \
                   self.sequence.count('C') + \
                   self.sequence.count('g') + \
                   self.sequence.count('c')
        return gc_count / len(self.sequence)

    def reverse_complement(self) -> str:
        """
        Returns the reverse complement of the sequence.

        Returns
        -------
        sequence   :   The reverse complement of the sequence.
        """
        complement = str.maketrans('ATCGatcg', 'TAGCtagc')
        return self.sequence.translate(complement)[::-1]

    def translate(self) -> PeptideSequence:
        """
        Translates the nucleotide sequence to a peptide sequence.

        Returns
        -------
        peptide_sequence:   PeptideSequence object.
        """
        nucleotide_seq = self.sequence.upper().replace('T', 'U')
        peptide_seq = ''
        for i in range(0, len(nucleotide_seq), 3):
            try:
                codon = nucleotide_seq[i:i + 3]
                amino_acid = AMINO_ACID_CODES[RNA_CODONS[codon]]
                peptide_seq += amino_acid
                aa_idx += 1
            except:
                raise Exception('')
        peptide_sequence = PeptideSequence(sequence=peptide_seq)
        return peptide_sequence
