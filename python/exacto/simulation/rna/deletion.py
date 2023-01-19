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
class Deletion(Exon):
    ref_exon: ClassVar[Exon]
    del_start: int = -1
    del_end: int = -1
    edits: List = field(default_factory=lambda: [])

    def __init__(self,
                 ref_exon: Exon,
                 del_start: int,
                 del_end: int):
        super().__init__()
        self.ref_exon = ref_exon
        self.del_start = del_start
        self.del_end = del_end
        self.edits = []

        # Append edits
        for i in range(0, len(self.ref_exon.edits)):
            if self.del_start <= self.ref_exon.edits[i].pos <= self.del_end:
                self.edits.append(Edit(
                    ref=self.ref_exon.edits[i].ref,
                    alt='',
                    pos=self.ref_exon.edits[i].pos,
                    sequence=''
                ))
            else:
                self.edits.append(self.ref_exon.edits[i])

    def __str__(self):
        msg = "[DELETION][%i:%i]\n" % (self.del_start, self.del_end)
        return msg + \
               super(Deletion, self).__str__()

