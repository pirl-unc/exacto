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
The purpose of this python3 script is to implement classes and functions
that are used to parse BAM files.
"""


import pysam
from typing import Tuple, List, Dict
from dataclasses import dataclass, field


def get_chrom_sizes(
        bam_file: pysam.AlignmentFile
    ) -> dict:
    """
    Returns chromosomes and their sizes in a BAM file.

    Parameters
    ----------
    bam_file            :   pysam.AlignmentFile object.

    Returns
    -------
    chrom_sizes_dict    :   Dictionary where keys are chromosomes
                            and values are their reference sizes.
    """
    chrom_sizes_dict = {} # key = chromosome, value = length
    for curr_chrom in bam_file.references:
        curr_size = bam_file.get_reference_length(curr_chrom)
        chrom_sizes_dict[curr_chrom] = curr_size
    return chrom_sizes_dict


def get_read_count(
        bam_file: pysam.AlignmentFile
    ) -> int:
    """
    Returns number of reads in a BAM file.

    Parameters
    ----------
    bam_file    :   pysam.AlignmentFile object.

    Returns
    -------
    num_reads   :   Number of reads
    """
    num_reads = 0
    for _ in bam_file.fetch():
        num_reads += 1
    return num_reads


def get_cs_tag(
        md_tag, cigar
    ) -> Tuple:
    """
    Returns the CS tag given a MD tag and a CIGAR string.

    Parameters
    ----------
    md_tag      :   MD tag.
    cigar       :   CIGAR string.

    Returns
    -------
    cs_tag      :
    """
    a = 1
