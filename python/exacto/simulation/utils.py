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


import random
import numpy as np
from typing import Tuple
from dataclasses import dataclass
from ..logging import get_logger
from ..constants import VariantTypes, Strands


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
    random_idx = random.randint(0, len(df_genes))
    df_curr_gene = df_genes.iloc[random_idx, :]
    curr_gene_id = df_curr_gene['gene_id'].values.tolist()[0]

    # Randomly select a transcript
    df_curr_gene_transcripts = df_transcripts.loc[df_transcripts['gene_id'] == curr_gene_id, :]
    random_idx = random.randint(0, len(df_curr_gene_transcripts))
    df_curr_transcript = df_curr_gene_transcripts.iloc[random_idx, :]
    curr_transcript_id = df_curr_transcript['transcript_id'].values.tolist()[0]
    curr_transcript_strand = df_curr_transcript['transcript_strand'].values.tolist()[0]

    # Randomly select an exon
    df_curr_transcript_exons = df_exons.loc[df_exons['transcript_id'] == curr_transcript_id, :]
    random_idx = random.randint(0, len(df_curr_transcript_exons))
    df_curr_transcript_exon = df_curr_transcript_exons.iloc[random_idx, :]

    # Randomly select a position
    pos = random.randint(df_curr_transcript_exon['exon_start'],
                         df_curr_transcript_exon['exon_end'])

    return RnaPos(
        gene_id=curr_gene_id,
        transcript_id=curr_transcript_id,
        exon_id=df_curr_transcript_exon['exon_id'].values.tolist()[0],
        exon_number=df_curr_transcript_exon['exon_number'].values.tolist()[0],
        exon_start=df_curr_transcript_exon['exon_start'].values.tolist()[0],
        exon_end=df_curr_transcript_exon['exon_end'].values.tolist()[0],
        pos=pos,
        chrom=df_curr_transcript_exon['exon_chrom'].values.tolist()[0],
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


def generate_insertion(insertion_size_mean: int,
                       insertion_size_stdev: int) -> str:
    """
    Generates an insertion.

    Parameters
    ----------
    insertion_size_mean     :   Mean value of insertion size.
    insertion_size_stdev    :   Standard deviation of insertion size.

    Returns
    -------
    insertion   :   A random nucleotide sequence.
    """
    size = int(np.random.normal(insertion_size_mean, insertion_size_stdev, 1)[0])
    atcg = ['A', 'T', 'C', 'G']
    sequence = [random.choice(atcg) for _ in range(0, size)]
    return ''.join(sequence)


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
    deletion_size = int(np.random.lognormal(mean=deletion_size_mean, sigma=deletion_size_stdev)[0])

    # Step 2. Manipulate the deletion position to fit the deletion type
    if deletion_type == VariantTypes.DeletionTypes.EXONIC_DELETION:
        if rna_pos.pos == rna_pos.exon_start:
            rna_pos.pos = rna_pos.pos + 1
        if rna_pos.pos == rna_pos.exon_end:
            rna_pos.pos = rna_pos.pos - 1
    if deletion_type == VariantTypes.DeletionTypes.FIVE_PRIME_SPLICE_SITE_DELETION:
        if rna_pos.strand == rna_pos.POSITIVE:
            rna_pos.pos = rna_pos.exon_start
        if rna_pos.strand == rna_pos.NEGATIVE:
            rna_pos.pos = rna_pos.exon_end
    if deletion_type == VariantTypes.DeletionTypes.THREE_PRIME_SPLICE_SITE_DELETION:
        if rna_pos.strand == Strands.POSITIVE:
            rna_pos.pos = rna_pos.exon_end
        if rna_pos.strand == Strands.NEGATIVE:
            rna_pos.pos = rna_pos.exon_start
    if rna_pos.strand == Strands.POSITIVE:
        if rna_pos.pos + deletion_size > rna_pos.exon_end:
            deletion_size = rna_pos.exon_end - rna_pos.pos + 1
    if rna_pos.strand == Strands.NEGATIVE:
        if rna_pos.pos - deletion_size < rna_pos.exon_start:
            deletion_size = rna_pos.pos - rna_pos.exon_start + 1

    # Step 3. Determine deletion start and end positions
    if rna_pos.strand == Strands.POSITIVE:
        deletion_start = rna_pos.pos
        deletion_end = rna_pos.pos + deletion_size - 1
    else:
        deletion_start = rna_pos.pos - deletion_size + 1
        deletion_end = rna_pos.pos

    return deletion_start, deletion_end, deletion_size


