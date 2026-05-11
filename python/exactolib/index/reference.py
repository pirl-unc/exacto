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
The purpose of this python3 script is to implement functions related to building indices for reference datasets.
"""


import pysam
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Set
from ..utilities import get_kmers


def _generate_kmers_for_chunk(args) -> Dict[int, Set[str]]:
    sequences, min_k, max_k = args
    result = {k: set() for k in range(min_k, max_k + 1)}
    for sequence in sequences:
        # Split on stop codons so no k-mer spans across one
        for fragment in sequence.split('*'):
            if not fragment:
                continue
            for k in range(min_k, max_k + 1):
                kmers = get_kmers(sequence=fragment, k=k)
                result[k].update(kmers)
    return result


def build_reference_kmer_index(
        reference_fasta_file: str,
        min_k: int,
        max_k: int,
        num_processes: int = 1
) -> Dict[int, Set[str]]:
    """
    Build a per-k set of all k-mers from a reference FASTA file.

    Args:
        reference_fasta_file    :   Path to reference FASTA file.
        min_k                   :   Minimum k-mer length.
        max_k                   :   Maximum k-mer length.
        num_processes           :   Number of processes.

    Returns:
        Dict mapping k (int) to Set[str] of all k-mers of that length.
    """
    # Step 1. Read all sequences into memory (sequential I/O)
    sequences: List[str] = []
    with pysam.FastxFile(reference_fasta_file) as fasta:
        for entry in fasta:
            sequences.append(str(entry.sequence))

    # Step 2. Generate k-mers (parallel compute)
    chunk_size = max(1, len(sequences) // num_processes)
    chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]

    with ProcessPoolExecutor(max_workers=num_processes) as pool:
        results = pool.map(_generate_kmers_for_chunk, [(chunk, min_k, max_k) for chunk in chunks])

    # Step 3. Merge results
    reference_kmers = {k: set() for k in range(min_k, max_k + 1)}
    for result in results:
        for k, kmers in result.items():
            reference_kmers[k].update(kmers)

    return reference_kmers
