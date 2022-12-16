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
from typing import Tuple, List


def generate_snv(reference_allele: str) -> str:
    """
    Generates a single-nucleotide variant.

    Parameters
    ----------
    reference_allele    :   Reference allele.

    Returns
    -------
    alternate_allele    :   Alternate allele.
    """
    reference_allele = reference_allele.upper()
    atcg = ['A', 'C', 'T', 'G']
    atcg.remove(reference_allele)
    alternate_allele = random.choice(atcg)
    return alternate_allele


def generate_random_nucleotide_sequence(size: int) -> str:
    """
    Generates a random nucleotide sequence.

    Parameters
    ----------
    size        :   Size.

    Returns
    -------
    insertion   :   A random nucleotide sequence.
    """
    atcg = ['A', 'T', 'C', 'G']
    sequence = [random.choice(atcg) for _ in range(0, size)]
    return ''.join(sequence)


def generate_deletion(size: int):
    """

    Parameters
    ----------
    size

    Returns
    -------

    """
    a = 1



