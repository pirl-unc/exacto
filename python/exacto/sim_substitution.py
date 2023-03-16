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


from dataclasses import dataclass
from typing import List, Tuple


from dataclasses import dataclass, field
from typing import List, Tuple, ClassVar
from .sim_edit import Edit
from .sim_exon import Exon


@dataclass
class Substitution(Exon):
    ref_exon: ClassVar[Exon]
    snv_pos: int = -1
    snv_alt: str = ''
    edits: List = field(default_factory=lambda: [])

    def __init__(
            self,
            ref_exon,
            snv_pos: int,
            snv_alt: str
        ):
        self.ref_exon = ref_exon
        self.snv_pos = snv_pos
        self.snv_alt = snv_alt
        self.edits = []

        # Append edits
        for i in range(0, len(self.ref_exon.edits)):
            if self.ref_exon.edits[i].pos == self.snv_pos:
                self.edits.append(Edit(
                    ref=self.ref_exon.edits[i].ref,
                    alt=self.snv_alt,
                    pos=self.ref_exon.edits[i].pos,
                    sequence=self.snv_alt
                ))
            else:
                self.edits.append(self.ref_exon.edits[i])

    def __str__(self):
        msg = "[SUBSTITUTION][%i:%i] %s\n" % (self.snv_pos, self.snv_pos, self.snv_alt)
        return msg + super(Substitution, self).__str__()
