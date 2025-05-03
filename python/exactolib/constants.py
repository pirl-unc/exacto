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
The purpose of this python3 script is to define Exacto constants.
"""


from enum import Enum, IntEnum


class DnaVariantCallingMode(Enum):
    CASE_SPECIFIC = 'case-specific'
    ALL = 'all'

    def __str__(self) -> str:
        return self.value


class GeneAnnotationSource(Enum):
    GENCODE = 'gencode'

    def __str__(self) -> str:
        return self.value


class OutputType(Enum):
    DATAFRAME = 'dataframe'
    FILE = 'file'

    def __str__(self) -> str:
        return self.value


class ReferenceTranscriptScoringMethod(Enum):
    COSINE_SIMILARITY = 'cosine_similarity'
    L2_DISTANCE = 'l2'
    JACCARD = 'jaccard'
    NET_OVERLAP = 'net_overlap'
    NON_OVERLAP = 'non_overlap'
    OVERLAP = 'overlap'
    WEIGHTED_NET_OVERLAP = 'weighted_net_overlap'

    def __str__(self) -> str:
        return self.value


class ReferenceTranscriptSelectionStrategy(Enum):
    TOP_K = 'top_k'
    THRESHOLD = 'threshold'

    def __str__(self) -> str:
        return self.value


class TranslationStrategy(Enum):
    LONGEST_ORF = 'longest-orf'

    def __str__(self) -> str:
        return self.value