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
The purpose of this python3 script is to implement the GraphGenomeMinimizer dataclasses.
"""


from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Set, Type


@dataclass(frozen=True)
class GraphGenomeBase:
    nucleotide: str
    node_id: int
    position: int


@dataclass(frozen=True)
class GraphGenomeSeqPos:
    node_id: int
    start: int
    end: int


@dataclass(frozen=True)
class GraphGenomeSeqAddr:
    positions: List[GraphGenomeSeqPos] = field(default_factory=list)

    @property
    def id(self):
        id = []
        for position in self.positions:
            id.append('%s:%i-%i' % (position.node_id, position.start, position.end))
        return ';'.join(id)


@dataclass(frozen=True)
class GraphGenomeMinimizerSet:
    k: int
    w: int
    minimizers: Dict[str, Dict[str, GraphGenomeSeqAddr]] = field(default_factory=lambda: defaultdict(dict))

    def add_minimizer(self, minimizer: str, seq_addr: GraphGenomeSeqAddr):
        if minimizer not in self.minimizers.keys():
            self.minimizers[minimizer] = {}

        if seq_addr.id not in self.minimizers[minimizer].keys():
            self.minimizers[minimizer][seq_addr.id] = seq_addr

