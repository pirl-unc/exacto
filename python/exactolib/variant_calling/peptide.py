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
The purpose of this python3 script is to implement classes and functions related to identifying peptide variants.
"""


import copy
import pandas as pd
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from typing import Dict, Set


@dataclass(frozen=False)
class AminoAcid:
    aa: str
    amino_acid_change: str
    primary_structure_index_start: int
    primary_structure_index_end: int
    rna_variant_call_ids: str
    dna_variant_call_ids: str


@dataclass(frozen=False)
class Peptide:
    peptide_id: int
    amino_acids: deque = field(default_factory=deque)

    def get_id(self) -> int:
        return self.peptide_id

    def get_length(self) -> int:
        return len(self.amino_acids)

    def get_sequence(self) -> str:
        return ''.join(aa.aa for aa in self.amino_acids)

    def push_back(self, amino_acid: AminoAcid):
        self.amino_acids.append(amino_acid)

    def pop_front(self):
        self.amino_acids.popleft()

    def is_mutant(self) -> bool:
        for aa in self.amino_acids:
            if aa.amino_acid_change == 'mutant':
                return True
        return False


def _worker(min_k, max_k, args):
    peptide_id, df_chunk = args

    # Keep base rows (remove event rows)
    df_chunk = df_chunk[df_chunk['type'] == 'base']

    mutant_peptides_dict = {} # key = length, value = List[Peptide]
    for k in range(min_k, max_k + 1):
        # Iterate through the peptide in the order of amino_acid index
        mutant_peptides = []
        peptide = Peptide(peptide_id=peptide_id)
        for _, group in df_chunk.groupby('amino_acid_index', sort=True):
            if len(group) != 3:
                break

            # Check if there is at least one mutant amino acid
            is_mutant = set(group['frameshift_state'].unique()) != {'inframe'}
            row0 = group.iloc[0]
            row1 = group.iloc[1]
            row2 = group.iloc[2]
            aa = str(row0['amino_acid'])

            if aa == '*':
                break

            primary_structure_index_start = row0['primary_structure_index']
            primary_structure_index_end = row2['primary_structure_index']

            rna_ids = set()
            dna_ids = set()
            for row in [row0, row1, row2]:
                val = row['rna_variant_call_ids']
                if isinstance(val, str) and val:
                    for id_str in val.split(','):
                        if id_str:
                            rna_ids.add(id_str)
                val = row['dna_variant_call_ids']
                if isinstance(val, str) and val:
                    for id_str in val.split(','):
                        if id_str:
                            dna_ids.add(id_str)
            rna_variant_call_ids = ','.join(sorted(rna_ids)) if rna_ids else ''
            dna_variant_call_ids = ','.join(sorted(dna_ids)) if dna_ids else ''

            assert row0['amino_acid_change'] == row1['amino_acid_change'] == row2['amino_acid_change'], 'Peptide ID: %s. Primary structure index: %i' % (peptide_id, primary_structure_index_start)

            amino_acid_change = row0['amino_acid_change']

            amino_acid = AminoAcid(
                aa=aa,
                amino_acid_change=amino_acid_change,
                primary_structure_index_start=primary_structure_index_start,
                primary_structure_index_end=primary_structure_index_end,
                rna_variant_call_ids=rna_variant_call_ids,
                dna_variant_call_ids=dna_variant_call_ids
            )

            peptide.push_back(amino_acid=amino_acid)

            if peptide.get_length() > k:
                peptide.pop_front()

            if peptide.get_length() == k and peptide.is_mutant():
                mutant_peptides.append(copy.deepcopy(peptide))

        mutant_peptides_dict[k] = mutant_peptides

    return mutant_peptides_dict


def _grouped_iter(df, col):
    for key in df[col].unique():
        yield key, df.loc[df[col] == key]


def identify_peptide_variants(
        df_primary_structures: pd.DataFrame,
        reference_kmer_set: Dict[int, Set[str]],
        min_k: int,
        max_k: int,
        num_processes: int
) -> pd.DataFrame:
    """
    Extract mutant peptide k-mers from primary structures and filter
    against a reference proteome.

    For each peptide_id, generates all k-mers (min_k to max_k) that overlap
    at least one mutant amino acid (frameshift_state != 'inframe'), then
    removes any k-mer found in the reference proteome.

    Args:
        df_primary_structures   :   DataFrame with the following columns:
        reference_kmer_set      :   Reference proteome k-mers (Dict[k, Set[k-mer]]).
        min_k                   :   Minimum peptide length.
        max_k                   :   Maximum peptide length.
        num_processes           :   Number of processes.

    Returns:
        pd.DataFrame with columns:
            peptide_id
            peptide_sequence
            primary_structure_index_start
            primary_structure_index_end,
            rna_variant_call_ids
            dna_variant_call_ids
    """
    # Step 1. Identify mutant peptides
    with ProcessPoolExecutor(max_workers=num_processes) as pool:
        results = list(pool.map(partial(_worker, min_k, max_k), _grouped_iter(df_primary_structures, "peptide_id")))

    # Step 2. Consolidate all the dictionaries into one
    mutant_peptides = defaultdict(list)
    for result in results:
        for k, peptides in result.items():
            mutant_peptides[k].extend(peptides)

    # Step 3. Identify mutant peptides absent in the reference peptides set
    sequence_to_id = {}
    next_id = 0
    data = {
        'mutant_peptide_id': [],
        'peptide_id': [],
        'mutant_peptide_sequence': [],
        'k': [],
        'primary_structure_index_start': [],
        'primary_structure_index_end': [],
        'rna_variant_call_ids': [],
        'dna_variant_call_ids': []
    }
    for k, peptides in mutant_peptides.items():
        ref_set = reference_kmer_set.get(k, set())
        for peptide in peptides:
            sequence = peptide.get_sequence()
            if sequence not in ref_set:
                if sequence not in sequence_to_id:
                    sequence_to_id[sequence] = next_id
                    next_id += 1
                rna_variant_call_ids = set()
                dna_variant_call_ids = set()
                for aa in peptide.amino_acids:
                    if isinstance(aa.rna_variant_call_ids, str) and aa.rna_variant_call_ids:
                        for id_str in aa.rna_variant_call_ids.split(','):
                            if id_str:
                                rna_variant_call_ids.add(id_str)
                    if isinstance(aa.dna_variant_call_ids, str) and aa.dna_variant_call_ids:
                        for id_str in aa.dna_variant_call_ids.split(','):
                            if id_str:
                                dna_variant_call_ids.add(id_str)
                data['peptide_id'].append(peptide.get_id())
                data['mutant_peptide_id'].append(sequence_to_id[sequence])
                data['mutant_peptide_sequence'].append(sequence)
                data['k'].append(len(sequence))
                data['primary_structure_index_start'].append(peptide.amino_acids[0].primary_structure_index_start)
                data['primary_structure_index_end'].append(peptide.amino_acids[-1].primary_structure_index_end)
                data['rna_variant_call_ids'].append(','.join(sorted(rna_variant_call_ids)))
                data['dna_variant_call_ids'].append(','.join(sorted(dna_variant_call_ids)))

    return pd.DataFrame(data)
