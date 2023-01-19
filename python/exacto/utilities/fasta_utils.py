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


from typing import Tuple, List, Dict
from dataclasses import dataclass, field


@dataclass
class Sequence:
    id: str
    sequence: str


def read_fasta_file(fasta_file: str) -> List[Sequence]:
    """
    Reads a FASTA file and returns sequence information.

    Parameters
    ----------
    fasta_file  :   FASTA file.

    Returns
    -------
    sequences   :   List of instances of the class Sequence.
    """
    sequences = []
    curr_id = ''
    with open(fasta_file, 'r') as f:
        for line in f.readlines():
            if '>' == line[0]:
                curr_id = line[1:].replace('\n', '')
            else:
                sequences.append(Sequence(id=curr_id, sequence=line.replace('\n','')))
    return sequences

