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


import random
import math
import gzip
import numpy as np
from collections import Counter
from typing import List
from dataclasses import dataclass
from exacto.fasta import Sequence
from exacto.logging import get_logger


logger = get_logger(__name__)


@dataclass
class Read:
    id: str
    sequence: str
    base_quality_score_string: str


def simulate_single_end_reads(
        sequences: List[Sequence],
        output_fastq_gz_file: str,
        num_bases: int,
        read_length_mean: float,
        read_length_stdev: float,
        base_quality_mean: float,
        base_quality_stdev: float
    ) -> None:
    """
    Simulates single-end reads.

    Parameters
    ----------
    sequences               :   List of instances of the class Sequence.
    output_fastq_gz_file    :   Output .fastq.gz file.
    num_bases               :   Number of bases to simulate.
    read_length_mean        :   Mean value of read length.
    read_length_stdev       :   Standard deviation of read length.
    base_quality_mean       :   Mean value of base quality.
    base_quality_stdev      :   Standard deviation of base quality.

    Returns
    -------
    reads               :   List of instance of the class Read.
    """
    # Step 1. Compute log-normal parameters based on normal distribution parameters
    read_length_lognormal_mean = math.log(math.pow(read_length_mean, 2) / math.sqrt(math.pow(read_length_mean, 2) + math.pow(read_length_stdev, 2)))
    read_length_lognormal_variance = math.log(1 + math.pow(read_length_stdev, 2) / math.pow(read_length_mean, 2))
    read_length_lognormal_sigma = math.sqrt(read_length_lognormal_variance)
    base_quality_lognormal_mean = math.log(math.pow(base_quality_mean, 2) / math.sqrt(math.pow(base_quality_mean, 2) + math.pow(base_quality_stdev, 2)))
    base_quality_lognormal_variance = math.log(1 + math.pow(base_quality_stdev, 2) / math.pow(base_quality_mean, 2))
    base_quality_lognormal_sigma = math.sqrt(base_quality_lognormal_variance)

    # Step 2. Simulate and write reads
    sequenced_molecules = Counter() # key = molecule ID, value = number of reads sequenced
    num_bases_simulated = 0
    num_reads_simulated = 0
    logger.info("Started writing simulated reads to FASTQ.GZ file.")
    with gzip.open(output_fastq_gz_file, 'wb') as f:
        while True:
            # Pick a sequence
            idx = random.randint(0, len(sequences) - 1)
            curr_molecule = sequences[idx]

            # Pick a read length
            curr_read_length = np.random.lognormal(mean=read_length_lognormal_mean, sigma=read_length_lognormal_sigma)

            # Read molecule sequence
            if curr_read_length >= len(curr_molecule.sequence):
                # Read the entire molecule sequence if read length is longer than molecule sequence length
                curr_sequence = curr_molecule.sequence
            else:
                # Read a part of the molecule sequence if read length is shorter than molecule sequence length
                idx = random.randint(0, len(curr_molecule.sequence) - curr_read_length)
                curr_sequence = curr_molecule.sequence[idx:]

            # Generate base quality scores
            curr_base_quality_scores = np.random.lognormal(mean=base_quality_lognormal_mean,
                                                           sigma=base_quality_lognormal_sigma,
                                                           size=len(curr_sequence))
            curr_base_quality_score = [chr(int(i)) for i in curr_base_quality_scores]
            curr_base_quality_score = ''.join(curr_base_quality_score)

            # Generate a read ID
            sequenced_molecules[curr_molecule.id] += 1
            read_id = '@' + str(curr_molecule.id) + '/' + str(sequenced_molecules[curr_molecule.id])

            # Write read to file
            f.write(str(read_id + '\n').encode())
            f.write(str(curr_sequence + '\n').encode())
            f.write(str('+' + '\n').encode())
            f.write(str(curr_base_quality_score + '\n').encode())
            num_reads_simulated += 1

            # Keep track of throughput
            num_bases_simulated += len(curr_sequence)
            if num_bases_simulated >= num_bases:
                break
    logger.info('%i reads were generated.' % num_reads_simulated)
    logger.info("Finished writing simulated reads to FASTQ.GZ file.")
