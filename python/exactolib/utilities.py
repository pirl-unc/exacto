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
The purpose of this python3 script is to implement general-purpose utility functions.
"""


import pysam
from typing import Dict, List, Set, Tuple
from .logging import get_logger


logger = get_logger(__name__)


def count_reads_in_fastx_file(fastx_file: str) -> int:
    count = 0
    with pysam.FastxFile(fastx_file) as fq:
        for _ in fq:
            count += 1
    return count


def get_chromosomes(bam_file: str) -> List[str]:
    """
    Get chromosomes in a BAM file.

    Parameters:
        bam_file    :   BAM file.

    Returns:
        List[chromosome]
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    return list(bam.references)


def get_chromosome_lengths(bam_file: str) -> Dict[str, int]:
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return dict(zip(bam.references, bam.lengths))


def get_chromosome_regions(bam_file: str) -> List[Tuple[str, int, int]]:
    """
    Get chromosome regions in the form <chrom>:1-<end> from a BAM file.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosomes = bam.references
        lengths = bam.lengths

    return [
        (chromosome, 1, length)   # tuple, not list
        for chromosome, length in zip(chromosomes, lengths)
    ]


def get_kmers(sequence: str, k: int) -> Set[str]:
    kmers = set()
    for i in range(0, len(sequence) - k + 1):
        kmers.add(sequence[i:i + k])
    return kmers


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise Exception('Boolean value expected.')

