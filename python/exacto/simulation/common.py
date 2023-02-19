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
The purpose of this python3 script is to implement common functions related
to the simulation module.
"""


import math
import numpy as np
import random
from typing import Tuple
from dataclasses import dataclass
from ..constants import VariantTypes, Strands
from ..logging import get_logger


logger = get_logger(__name__)


@dataclass
class RnaPos:
    gene_id: str = ''
    transcript_id: str = ''
    exon_id: str = ''
    exon_number: int = -1
    exon_start: int = -1
    exon_end: int = -1
    pos: int = -1
    chrom: str = ''
    strand: str = ''


def randomly_select_rna_position(df_genes,
                                 df_transcripts,
                                 df_exons) -> RnaPos:
    """
    Randomly selects a RNA position.

    Parameters
    ----------
    df_genes            :   DataFrame of reference genes.
    df_transcripts      :   DataFrame of reference transcripts.
    df_exons            :   DataFrame of reference exons.

    Returns
    -------
    random_pos          :   an instance of RandomPos
    """
    # Randomly select a gene
    random_idx = random.randint(0, len(df_genes) - 1)
    df_curr_gene = df_genes.iloc[random_idx, :]
    curr_gene_id = df_curr_gene['gene_id']

    # Randomly select a transcript
    df_curr_gene_transcripts = df_transcripts.loc[df_transcripts['gene_id'] == curr_gene_id,:]
    random_idx = random.randint(0, len(df_curr_gene_transcripts) - 1)
    df_curr_transcript = df_curr_gene_transcripts.iloc[random_idx, :]
    curr_transcript_id = df_curr_transcript['transcript_id']
    curr_transcript_strand = df_curr_transcript['transcript_strand']

    # Randomly select an exon
    df_curr_transcript_exons = df_exons.loc[df_exons['transcript_id'] == curr_transcript_id,:]
    random_idx = random.randint(0, len(df_curr_transcript_exons) - 1)
    df_curr_transcript_exon = df_curr_transcript_exons.iloc[random_idx,:]

    # Randomly select a position
    pos = random.randint(df_curr_transcript_exon['exon_start'],
                         df_curr_transcript_exon['exon_end'])

    return RnaPos(
        gene_id=curr_gene_id,
        transcript_id=curr_transcript_id,
        exon_id=df_curr_transcript_exon['exon_id'],
        exon_number=df_curr_transcript_exon['exon_number'],
        exon_start=df_curr_transcript_exon['exon_start'],
        exon_end=df_curr_transcript_exon['exon_end'],
        pos=pos,
        chrom=df_curr_transcript_exon['exon_chrom'],
        strand=curr_transcript_strand
    )


def generate_single_nucleotide_variant(reference_allele: str) -> str:
    """
    Generates a single-nucleotide RNA variant.

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


def generate_random_sequence(size, disallowed_sequences):
    """
    Generate a random sequence.

    Parameters
    ----------
    size                    :   Size of random sequence.
    disallowed_sequences    :   List of disallowed sequences.

    Returns
    -------
    random_sequence         :   Random sequence.
    """
    atcg = ['A', 'T', 'C', 'G']
    while True:
        random_sequence = [random.choice(atcg) for _ in range(0, size)]
        random_sequence = ''.join(random_sequence)
        allowed = True
        for seq in disallowed_sequences:
            if random_sequence in seq:
                allowed = False
        if allowed:
            break
    return random_sequence


def generate_insertion(insertion_size_mean: int,
                       insertion_size_stdev: int,
                       disallowed_sequences: []) -> str:
    """
    Generates an insertion.

    Parameters
    ----------
    insertion_size_mean      :   Mean value of insertion size.
    insertion_size_stdev     :   Standard deviation of insertion size.
    disallowed_sequences     :   List of disallowed subsequences.

    Returns
    -------
    insertion   :   A random nucleotide sequence.
    """
    lognormal_mean = math.log(math.pow(insertion_size_mean, 2) / math.sqrt(math.pow(insertion_size_mean, 2) + math.pow(insertion_size_stdev, 2)))
    lognormal_variance = math.log(1 + math.pow(insertion_size_stdev, 2) / math.pow(insertion_size_mean, 2))
    lognormal_sigma = math.sqrt(lognormal_variance)
    insertion_size = int(np.random.lognormal(
        mean=lognormal_mean,
        sigma=lognormal_sigma
    ))
    atcg = ['A', 'T', 'C', 'G']
    insertion_sequence = [random.choice(atcg) for _ in range(0, insertion_size)]
    insertion_sequence = ''.join(insertion_sequence)
    disallowed_subsequences = [i.upper() for i in disallowed_sequences]

    # Find occurrences of all disallowed sequences and change them
    while True:
        # For each disallowed sequence, replace it with an allowed sequence
        for subseq in disallowed_subsequences:
            subseq_positions = [i for i in range(len(insertion_sequence)) if insertion_sequence.startswith(subseq, i)]
            for curr_pos in subseq_positions:
                # Generate a new sequence
                new_subsequence = generate_random_sequence(size=len(subseq), disallowed_sequences=disallowed_sequences)
                insertion_sequence[:curr_pos] + new_subsequence + insertion_sequence[curr_pos + len(subseq):]

        # Make sure none of the disallowed sequences exists
        allowed = True
        for subseq in disallowed_subsequences:
            subseq_positions = [i for i in range(len(insertion_sequence)) if insertion_sequence.startswith(subseq, i)]
            if len(subseq_positions) > 0:
                allowed = False
                break
        if allowed:
            break

    return insertion_sequence


def generate_rna_deletion(rna_pos: RnaPos,
                          deletion_size_mean: int,
                          deletion_size_stdev: int) -> Tuple[int, int, int]:
    """
    Generates a deletion. Deletion size is sampled from a log-normal distribution.

    Parameters
    ----------
    rna_pos                 :   An instance of the class RnaPos.
    deletion_size_mean      :   Mean value of deletion size.
    deletion_size_stdev     :   Standard deviation of deletion size.

    Returns
    -------
    deletion_start          :   Deletion start position.
    deletion_end            :   Deletion end position.
    deletion_size           :   Deletion size.
    """
    # Step 1. Random select a deletion type
    deletion_type = random.choice(VariantTypes.DeletionTypes.ALL)
    lognormal_mean = math.log(math.pow(deletion_size_mean, 2) / math.sqrt(math.pow(deletion_size_mean, 2) + math.pow(deletion_size_stdev, 2)))
    lognormal_variance = math.log(1 + math.pow(deletion_size_stdev, 2) / math.pow(deletion_size_mean, 2))
    lognormal_sigma = math.sqrt(lognormal_variance)
    deletion_size = int(np.random.lognormal(
        mean=lognormal_mean,
        sigma=lognormal_sigma
    ))

    # Step 2. Manipulate the deletion position to fit the deletion type
    if deletion_type == VariantTypes.DeletionTypes.EXONIC_DELETION:
        if rna_pos.pos == rna_pos.exon_start:
            rna_pos.pos = rna_pos.pos + 1
        if rna_pos.pos == rna_pos.exon_end:
            rna_pos.pos = rna_pos.pos - 1
    if deletion_type == VariantTypes.DeletionTypes.FIVE_PRIME_SPLICE_SITE_DELETION:
        if rna_pos.strand == Strands.POSITIVE:
            rna_pos.pos = rna_pos.exon_start
        if rna_pos.strand == Strands.NEGATIVE:
            rna_pos.pos = rna_pos.exon_end
    if deletion_type == VariantTypes.DeletionTypes.THREE_PRIME_SPLICE_SITE_DELETION:
        if rna_pos.strand == Strands.POSITIVE:
            rna_pos.pos = rna_pos.exon_end
        if rna_pos.strand == Strands.NEGATIVE:
            rna_pos.pos = rna_pos.exon_start
    if rna_pos.strand == Strands.POSITIVE:
        if (rna_pos.pos + deletion_size - 1) > rna_pos.exon_end:
            deletion_size = rna_pos.exon_end - rna_pos.pos + 1
    if rna_pos.strand == Strands.NEGATIVE:
        if (rna_pos.pos - deletion_size + 1) < rna_pos.exon_start:
            deletion_size = rna_pos.pos - rna_pos.exon_start + 1

    # Step 3. Determine deletion start and end positions
    if rna_pos.strand == Strands.POSITIVE:
        deletion_start = rna_pos.pos
        deletion_end = rna_pos.pos + deletion_size - 1
    else:
        deletion_start = rna_pos.pos - deletion_size + 1
        deletion_end = rna_pos.pos

    return deletion_start, deletion_end, deletion_size


