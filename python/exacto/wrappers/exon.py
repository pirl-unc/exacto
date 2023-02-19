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


class Exon:

    def __init__(
            self,
            ref_id: str,
            ref_start: int,
            ref_end: int,
            ref_sequence: str
        ):
        self.__ref_id = ref_id
        self.__ref_start = ref_start
        self.__ref_end = ref_end
        self.__ref_sequence = ref_sequence
        self.__sequence = "" # actual exon sequence

        # Linked list variables
        self.__next_exon = None

    def get_ref_id(self):
        return self.__ref_id

    def get_ref_start(self):
        return self.__ref_start

    def get_ref_end(self):
        return self.__ref_end

    def get_ref_sequence(self):
        return self.__ref_sequence

    def get_sequence(self):
        return self.__sequence

    def set_sequence(self, sequence):
        self.__sequence = sequence

    def get_next_exon(self):
        return self.__next_exon

    def set_next_exon(self, exon):
        self.__next_exon = exon
