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
Identify peptide variants from primary structures.

The input DataFrame has one row per primary structure (`PrimaryStructureRecord`):
its `amino_acid_sequence` column carries the full peptide as a string,
`mutant_amino_acid_intervals` encodes the mutant positions as a compact
interval string (e.g. `"3-5,12,20-21"`, 0-indexed inclusive), and
`rna_variant_ids` / `dna_variant_ids` list the contributing variants at the
primary-structure level.

For each primary structure we slide a window of size k over the amino-acid
sequence, emit every window that overlaps at least one mutant position
(skipping windows that contain a stop codon `*`), and drop any k-mer found
in the reference proteome.
"""


import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Dict, List, Set, Tuple


# Per-row payload threaded through the worker pool. Kept as a tuple (not a
# dataclass) so it pickles cheaply for multiprocessing.
PrimaryStructureRow = Tuple[int, str, str, str, str]

# Worker output: (primary_structure_id, kmer, aa_index_start, aa_index_end,
#                 rna_variant_ids, dna_variant_ids).
MutantKmer = Tuple[int, str, int, int, str, str]


def _parse_intervals(intervals_str: str) -> Set[int]:
    """Parse a compact interval string into the set of integers it covers.

    Accepts the same format `format_intervals` emits on the Rust side:
    comma-separated runs that are either a single number `X` or an inclusive
    range `X-Y`. Empty / non-string input yields an empty set.
    """
    if not isinstance(intervals_str, str) or not intervals_str:
        return set()
    positions: Set[int] = set()
    for part in intervals_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            start_str, end_str = part.split('-', 1)
            positions.update(range(int(start_str), int(end_str) + 1))
        else:
            positions.add(int(part))
    return positions


def _worker(min_k: int, max_k: int, row: PrimaryStructureRow) -> List[MutantKmer]:
    """Extract every mutant k-mer (min_k <= k <= max_k) from one primary structure.

    A k-mer is "mutant" iff its window overlaps at least one mutant amino-acid
    index. K-mers containing the stop codon `*` are skipped.
    """
    primary_structure_id, sequence, intervals_str, rna_variant_ids, dna_variant_ids = row

    if not isinstance(sequence, str) or not sequence:
        return []

    mutant_positions = _parse_intervals(intervals_str)
    if not mutant_positions:
        return []

    rna_variant_ids = rna_variant_ids if isinstance(rna_variant_ids, str) else ''
    dna_variant_ids = dna_variant_ids if isinstance(dna_variant_ids, str) else ''

    results: List[MutantKmer] = []
    seq_len = len(sequence)
    for k in range(min_k, max_k + 1):
        if k > seq_len:
            break
        for start in range(seq_len - k + 1):
            end = start + k - 1
            # Cheap overlap check: do any mutant positions land in [start, end]?
            if not any((start <= p <= end) for p in mutant_positions):
                continue
            kmer = sequence[start:start + k]
            if '*' in kmer:
                continue
            results.append((primary_structure_id, kmer, start, end, rna_variant_ids, dna_variant_ids))
    return results


def identify_peptide_variants(
        df_primary_structures: pd.DataFrame,
        reference_kmer_set: Dict[int, Set[str]],
        min_k: int,
        max_k: int,
        num_processes: int
) -> pd.DataFrame:
    """
    Extract mutant peptide k-mers from primary structures and filter them
    against a reference proteome.

    Args:
        df_primary_structures   :   DataFrame with one row per primary
                                    structure. Required columns:
                                        primary_structure_id,
                                        amino_acid_sequence,
                                        mutant_amino_acid_intervals,
                                        rna_variant_ids,
                                        dna_variant_ids.
        reference_kmer_set      :   Reference proteome k-mers
                                    (Dict[k, Set[k-mer]]).
        min_k                   :   Minimum peptide length.
        max_k                   :   Maximum peptide length.
        num_processes           :   Number of worker processes.

    Returns:
        pd.DataFrame with columns:
            mutant_peptide_id,
            primary_structure_id,
            mutant_peptide_sequence,
            k,
            amino_acid_index_start,
            amino_acid_index_end,
            rna_variant_ids,
            dna_variant_ids
    """
    # Step 1. Assemble per-row payloads. zip over Series is faster than
    # iterrows() and avoids the per-row dtype-coerce cost.
    row_payloads: List[PrimaryStructureRow] = list(zip(
        df_primary_structures['primary_structure_id'].astype(int).tolist(),
        df_primary_structures['amino_acid_sequence'].fillna('').astype(str).tolist(),
        df_primary_structures['mutant_amino_acid_intervals'].fillna('').astype(str).tolist(),
        df_primary_structures['rna_variant_ids'].fillna('').astype(str).tolist(),
        df_primary_structures['dna_variant_ids'].fillna('').astype(str).tolist(),
    ))

    # Step 2. Extract mutant k-mers in parallel.
    if num_processes <= 1 or len(row_payloads) <= 1:
        nested_results: List[List[MutantKmer]] = [
            _worker(min_k, max_k, row) for row in row_payloads
        ]
    else:
        with ProcessPoolExecutor(max_workers=num_processes) as pool:
            nested_results = list(pool.map(partial(_worker, min_k, max_k), row_payloads))

    # Step 3. Flatten, filter against the reference proteome, dedupe sequences.
    sequence_to_id: Dict[str, int] = {}
    next_id = 0
    data: Dict[str, list] = {
        'mutant_peptide_id': [],
        'primary_structure_id': [],
        'mutant_peptide_sequence': [],
        'k': [],
        'amino_acid_index_start': [],
        'amino_acid_index_end': [],
        'rna_variant_ids': [],
        'dna_variant_ids': []
    }
    for results in nested_results:
        for (primary_structure_id, sequence, aa_start, aa_end, rna_ids, dna_ids) in results:
            k = len(sequence)
            ref_set = reference_kmer_set.get(k, set())
            if sequence in ref_set:
                continue
            if sequence not in sequence_to_id:
                sequence_to_id[sequence] = next_id
                next_id += 1
            data['mutant_peptide_id'].append(sequence_to_id[sequence])
            data['primary_structure_id'].append(primary_structure_id)
            data['mutant_peptide_sequence'].append(sequence)
            data['k'].append(k)
            data['amino_acid_index_start'].append(aa_start)
            data['amino_acid_index_end'].append(aa_end)
            data['rna_variant_ids'].append(rna_ids)
            data['dna_variant_ids'].append(dna_ids)

    return pd.DataFrame(data)
