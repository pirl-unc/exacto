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


import pandas as pd
from .exon import Exon


class Transcript:

    def __init__(self,
                 ref_id: str,
                 ref_chromosome: str,
                 ref_start: int,
                 ref_end: int,
                 ref_type: str,
                 ref_strand: str):
        # Reference (original) transcript information
        self.__ref_id = ref_id
        self.__ref_chromosome = ref_chromosome
        self.__ref_start = ref_start
        self.__ref_end = ref_end
        self.__ref_type = ref_type
        self.__ref_strand = ref_strand

        # Linked list head exon
        self.__head_exon = None

    def get_ref_id(self):
        return self.__ref_id

    def get_ref_chromosome(self):
        return self.__ref_chromosome

    def get_ref_start(self):
        return self.__ref_start

    def get_ref_end(self):
        return self.__ref_end

    def get_ref_type(self):
        return self.__ref_type

    def get_ref_strand(self):
        return self.__ref_strand

    def get_head_exon(self):
        return self.__head_exon

    def insert_exon(self, exon):
        if self.__head_exon:
            current = self.__head_exon
            while current.get_next_exon():
                current = current.get_next_exon()
            current.set_next_exon(exon=exon)
        else:
            self.__head_exon = exon

