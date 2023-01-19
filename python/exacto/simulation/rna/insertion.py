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
from typing import List, Tuple, ClassVar
from .edit import Edit
from .exon import Exon


@dataclass
class Insertion(Exon):
    ref_exon: ClassVar[Exon]
    ins_pos: int = -1
    ins_sequence: str = ''
    edits: List = field(default_factory=lambda: [])

    def __init__(self,
                 ref_exon,
                 ins_pos: int,
                 ins_sequence: str):
        self.ref_exon = ref_exon
        self.ins_pos = ins_pos
        self.ins_sequence = ins_sequence
        self.edits = []

        # Append edits
        for i in range(0, len(self.ref_exon.edits)):
            if self.ref_exon.edits[i].pos == self.ins_pos:
                self.edits.append(self.ref_exon.edits[i])
                self.edits.append(Edit(
                    ref='',
                    alt=self.ins_sequence,
                    pos=self.ref_exon.edits[i].pos,
                    sequence=self.ins_sequence
                ))
            else:
                self.edits.append(self.ref_exon.edits[i])

    def __str__(self):
        msg = "[INSERTION][%i:%i] %s\n" % (self.ins_pos, self.ins_pos, self.ins_sequence)
        return msg + \
               super(Insertion, self).__str__()
