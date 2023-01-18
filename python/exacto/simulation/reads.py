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
import time
import numpy as np
from collections import Counter
from typing import Tuple, List, Dict
from dataclasses import dataclass, field
from ..utilities.fasta_utils import Sequence
from ..logging import get_logger


logger = get_logger(__name__)


@dataclass
class Read:
    id: str
    sequence: str
    base_quality_score_string: str


def simulate_single_end_reads(sequences: List[Sequence],
                              num_bases: int,
                              read_length_mean: float,
                              read_length_stdev: float,
                              base_quality_mean: float,
                              base_quality_stdev: float) -> List[Read]:
    """
    Simulates single-end reads.

    Parameters
    ----------
    sequences           :   List of instances of the class Sequence.
    num_bases           :   Number of bases to simulate.
    read_length_mean    :   Mean value of read length.
    read_length_stdev   :   Standard deviation of read length.
    base_quality_mean   :   Mean value of base quality.
    base_quality_stdev  :   Standard deviation of base quality.

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

    # Step 2. Simulate reads
    reads = []
    sequenced_molecules = Counter() # key = molecule ID, value = number of reads sequenced
    num_bases_simulated = 0
    while True:
        # Pick a sequence
        idx = random.randint(0, len(sequences) - 1)
        curr_molecule = sequences[idx]

        # Pick a read length
        curr_read_length = np.random.lognormal(mean=read_length_lognormal_mean, sigma=read_length_lognormal_sigma)

        # Read molecule sequence
        if curr_read_length >= len(curr_molecule.sequence):
            curr_sequence = curr_molecule.sequence
        else:
            # Pick a location of the sequence
            idx = random.randint(0, len(curr_molecule.sequence) - curr_read_length)
            curr_sequence = curr_molecule.sequence[idx:]

        # Generate base quality scores
        curr_base_quality_score = ''
        for i in range(0, len(curr_sequence)):
            curr_score = np.random.lognormal(mean=base_quality_lognormal_mean, sigma=base_quality_lognormal_sigma)
            curr_base_quality_score = curr_base_quality_score + chr(int(curr_score))

        # Generate a read ID
        sequenced_molecules[curr_molecule.id] += 1
        read_id = '@' + str(curr_molecule.id) + '/' + str(sequenced_molecules[curr_molecule.id])

        # Store sequenced read
        read = Read(id=read_id,
                    sequence=curr_sequence,
                    base_quality_score_string=curr_base_quality_score)
        reads.append(read)

        num_bases_simulated += len(curr_sequence)
        if num_bases_simulated >= num_bases:
            break

    # Provide information on how many molecules were sequenced
    logger.info('At least 1 read was generated for %i/%i molecules (sequences).'
                % (len(sequenced_molecules.keys()), len(sequences)))

    return reads

