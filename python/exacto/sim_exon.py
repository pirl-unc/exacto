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


@dataclass
class Exon:
    edits: list

    def __init__(self):
        self.edits = []

    @property
    def sequence(self):
        return ''.join([edit.sequence for edit in self.edits])

    def __str__(self):
        values = ['\tedits\t\t\t:']
        for i in range(0, len(self.edits)):
            if i == 0:
                values.append('\t')
            else:
                values.append('\t\t\t\t\t')
            if self.edits[i].type == 'REF':
                values.append('[%s][%i:%s]' %
                              (self.edits[i].type,
                               self.edits[i].pos,
                               self.edits[i].ref))
            if self.edits[i].type == 'SNV':
                values.append('[%s][%i:%s>%s]' %
                              (self.edits[i].type,
                               self.edits[i].pos,
                               self.edits[i].ref,
                               self.edits[i].alt))
            if self.edits[i].type == 'DEL':
                values.append('[%s][%i:%s]' %
                              (self.edits[i].type,
                               self.edits[i].pos,
                               self.edits[i].ref))
            if self.edits[i].type == 'INS':
                values.append('[%s][%i:%s]' %
                              (self.edits[i].type,
                               self.edits[i].pos,
                               self.edits[i].alt))
            values.append('\n')
        values.append('\tcurrent sequence\t:\t%s' % self.sequence)
        return ''.join(values)
