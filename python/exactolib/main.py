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
The purpose of this python3 script is to implement Exacto's main APIs.
"""


import json
import pysam
from typing import List, Tuple
from .default import *
from .logging import get_logger
from .utilities import get_chromosomes
from .variant_call import VariantCall
from exactolib import exactolibrs


logger = get_logger(__name__)


def identify_rna_variants(
        bam_file: str,
        min_reads: int,
        min_mapping_quality: int,
        num_threads: int,
        min_ins_size_proportion: float = IDENTIFY_MIN_INS_SIZE_PROPORTION,
        max_ins_norm_edit_distance: float = IDENTIFY_MAX_INS_NORM_EDIT_DISTANCE,
        min_del_size_proportion: float = IDENTIFY_MIN_DEL_SIZE_PROPORTION,
        chromosomes: List[str] = None
) -> List[VariantCall]:
    """
    Call RNA variants in a long-read RNA-seq BAM file.

    Parameters:
        bam_file                    :   BAM file.
        min_reads                   :   Minimum number of reads.
        min_mapping_quality         :   Minimum mapping quality.
        num_threads                 :   Number of threads.
        min_ins_size_proportion     :   Minimum insertion size proportion between
                                        two insertions. Size proportion = smaller
                                        insertion size / longer insertion size.
        max_ins_norm_edit_distance  :   Maximum insertion normalized edit
                                        (Levenshtein) distance. Normalized edit
                                        distance = edit distance / longer insertion size.
        min_del_size_proportion     :   Minimum deletion size proportion between
                                        two deletions. Size proportion = smaller
                                        deletion size / longer deletion size.
        chromosomes                 :   Chromosomes to call variants
                                        (if unspecified, variants are called in all
                                        chromosomes).

    Returns:
        List[VariantCall]
    """
    if chromosomes is None:
        chromosomes = get_chromosomes(bam_file=bam_file)
    json_str = exactolibrs.identify_rna_variants(
        bam_file,
        min_reads,
        min_mapping_quality,
        num_threads,
        min_ins_size_proportion,
        max_ins_norm_edit_distance,
        min_del_size_proportion,
        chromosomes
    )
    variant_calls = []
    for data in json.loads(json_str):
        variant_call = VariantCall(**data)
        variant_calls.append(variant_call)
    return variant_calls
