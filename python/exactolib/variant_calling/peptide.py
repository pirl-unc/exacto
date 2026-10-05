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


import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Dict, List, Set, Tuple


# Per-row payload threaded through the worker pool. Kept as a tuple (not a
# dataclass) so it pickles cheaply for multiprocessing.
ProteoformRow = Tuple[int, str, str, str, str]

# Worker output: (proteoform_id, kmer, aa_index_start, aa_index_end,
#                 assembled_transcript_variant_ids, dna_variant_ids).
MutantKmer = Tuple[int, str, int, int, str, str]


def _parse_intervals(intervals_str: str) -> Set[int]:
    if not isinstance(intervals_str, str) or not intervals_str:
        return set()
    positions: Set[int] = set()
    for part in re.split(r'[;,]', intervals_str):
        part = part.strip()
        if not part:
            continue
        bounds = re.split(r'[:-]', part)
        if len(bounds) == 2:
            positions.update(range(int(bounds[0]), int(bounds[1]) + 1))
        else:
            positions.add(int(part))
    return positions


def _worker(min_k: int, max_k: int, row: ProteoformRow) -> List[MutantKmer]:
    proteoform_id, sequence, intervals_str, assembled_transcript_variant_ids, dna_variant_ids = row

    if not isinstance(sequence, str) or not sequence:
        return []

    mutant_positions = _parse_intervals(intervals_str)
    if not mutant_positions:
        return []

    assembled_transcript_variant_ids = assembled_transcript_variant_ids if isinstance(assembled_transcript_variant_ids, str) else ''
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
            results.append((
                proteoform_id, kmer, start, end,
                assembled_transcript_variant_ids, dna_variant_ids
            ))
    return results


def identify_peptide_variants(
        df_proteoforms: pd.DataFrame,
        reference_kmer_set: Dict[int, Set[str]],
        min_k: int,
        max_k: int,
        num_processes: int
) -> pd.DataFrame:
    """
    Extract mutant peptide k-mers from proteoforms and filter them against a
    reference proteome.

    Args:
        df_proteoforms          :   DataFrame with one row per proteoform, as
                                    written by translate-transcripts. Required
                                    columns:
                                        proteoform_id,
                                        amino_acid_sequence,
                                        mutant_amino_acid_intervals,
                                        assembled_transcript_variant_ids,
                                        dna_variant_ids.
        reference_kmer_set      :   Reference proteome k-mers
                                    (Dict[k, Set[k-mer]]).
        min_k                   :   Minimum peptide length.
        max_k                   :   Maximum peptide length.
        num_processes           :   Number of worker processes.

    Returns:
        pd.DataFrame with columns:
            mutant_peptide_id,
            proteoform_id,
            mutant_peptide_sequence,
            k,
            amino_acid_index_start,
            amino_acid_index_end,
            assembled_transcript_variant_ids,
            dna_variant_ids
    """
    # Step 1. Assemble per-row payloads.
    row_payloads: List[ProteoformRow] = list(zip(
        df_proteoforms['proteoform_id'].astype(int).tolist(),
        df_proteoforms['amino_acid_sequence'].fillna('').astype(str).tolist(),
        df_proteoforms['mutant_amino_acid_intervals'].fillna('').astype(str).tolist(),
        df_proteoforms['assembled_transcript_variant_ids'].fillna('').astype(str).tolist(),
        df_proteoforms['dna_variant_ids'].fillna('').astype(str).tolist(),
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
        'proteoform_id': [],
        'mutant_peptide_sequence': [],
        'k': [],
        'amino_acid_index_start': [],
        'amino_acid_index_end': [],
        'assembled_transcript_variant_ids': [],
        'dna_variant_ids': []
    }
    for results in nested_results:
        for (proteoform_id, sequence, aa_start, aa_end, assembled_transcript_ids, dna_ids) in results:
            k = len(sequence)
            ref_set = reference_kmer_set.get(k, set())
            if sequence in ref_set:
                continue
            if sequence not in sequence_to_id:
                sequence_to_id[sequence] = next_id
                next_id += 1
            data['mutant_peptide_id'].append(sequence_to_id[sequence])
            data['proteoform_id'].append(proteoform_id)
            data['mutant_peptide_sequence'].append(sequence)
            data['k'].append(k)
            data['amino_acid_index_start'].append(aa_start)
            data['amino_acid_index_end'].append(aa_end)
            data['assembled_transcript_variant_ids'].append(assembled_transcript_ids)
            data['dna_variant_ids'].append(dna_ids)

    return pd.DataFrame(data)
