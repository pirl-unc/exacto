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


from typing import Tuple, List
from .transcript import Transcript


class Gene:

    def __init__(
            self,
            ref_id: str,
            ref_chromosome: str,
            ref_start: int,
            ref_end: int,
            ref_strand: str
        ):
        self.__ref_id = ref_id
        self.__ref_chromosome = ref_chromosome
        self.__ref_start = ref_start
        self.__ref_end = ref_end
        self.__ref_strand = ref_strand # '+' or '-'
        self.__transcripts = [] # list of Transcript objects

    def get_ref_id(self) -> str:
        return self.__ref_id

    def get_ref_chromosome(self) -> str:
        return self.__ref_chromosome

    def get_ref_start(self) -> int:
        return self.__ref_start

    def get_ref_end(self) -> int:
        return self.__ref_end

    def get_ref_strand(self) -> str:
        return self.__ref_strand

    def get_transcripts(self) -> List[Transcript]:
        return self.__transcripts

    def add_transcript(
            self,
            transcript
        ):
        """
        Adds a transcript.

        Args
        ----
        transcript  :   Transcript object.
        """
        self.__transcripts.append(transcript)
