#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
simulating transcript-level single-nucleotide variants,
small insertions and deletions.

Author: Jin Seok (Andy) Lee

Last updated date: July 20, 2022
"""


import random
import pandas as pd
import numpy as np
import logging
from typing import Tuple
from pysam import FastaFile
from exactolib.simulation.common import *


def randomly_generate_snv(transcript_id: str,
                          df_gtf: pd.DataFrame,
                          reference_genome_fasta_file: str) -> Tuple[str, str, str, str, int, int, str, str]:
    """
    Randomly generates a single-nucleotide variants.

    Args
    ----
    transcript_id               :   Transcript ID.
    df_gtf                      :   This DataFrame should only contain all the exons
                                    of one transcript. DataFrame with the following columns:
                                    'transcript_id',
                                    'exon_id',
                                    'exon_number',
                                    'chrom',
                                    'start' (reference start position),
                                    'end' (reference end position),
                                    'size'.
    reference_genome_fasta_file :   Reference genome FASTA file.

    Returns
    -------
    variant_type,
    strand,
    variant_transcript_sequence,
    variant_chrom,
    variant_position (reference position),
    exon_number,
    ref_allele,
    var_allele
    """
    # Step 1. Get exons that belong to transcript_id
    df_gtf_transcript = df_gtf.loc[df_gtf['transcript_id'] == transcript_id, :]

    # Step 2. Sort by exon_number
    df_gtf_transcript = df_gtf_transcript.sort_values(by='exon_number', ascending=True)

    # Step 3. Randomly select an exon
    exon_number = random.choice(df_gtf_transcript['exon_number'].values.tolist())

    # Step 4. Iterate through each exon and generate a variant transcript sequence
    fasta_object = FastaFile(reference_genome_fasta_file)
    variant_transcript_sequence = ''
    for index, row in df_gtf_transcript.iterrows():
        sequence = fasta_object.fetch(row['chrom'], row['start'] - 1, row['end'])
        if row['exon_number'] == exon_number:
            # Randomly select a position on the exon
            variant_position_local = random.randint(0, len(sequence))
            variant_position_ref = row['start'] + variant_position_local

            # Randomly select a variant allele
            ref_allele = fasta_object.fetch(row['chrom'],
                                            variant_position_ref - 1,
                                            variant_position_ref)
            atcg = ['A', 'T', 'C', 'G']
            atcg.remove(ref_allele)
            var_allele = random.choice(atcg)

            # Piece together the variant sequence
            sequence = sequence[0:variant_position_local] + \
                       var_allele + \
                       sequence[variant_position_local + 1:]
        variant_transcript_sequence += sequence

    # Step 5. Reverse the sequence if on the minus strand
    strand = df_gtf_transcript['strand'].unique()[0]
    if strand == '-':
        variant_transcript_sequence = variant_transcript_sequence[::-1]
    variant_chrom = df_gtf_transcript['chrom'].unique()[0]

    return 'snv', \
           strand, \
           variant_transcript_sequence, \
           variant_chrom, \
           variant_position_ref, \
           exon_number, \
           ref_allele, \
           var_allele


def randomly_generate_indel(sequence: str,
                            strand: str,
                            df_exons: pd.DataFrame) -> Tuple[str, str, int, int]:
    """
    Randomly generates a small insertion or deletion (< 50-bp).

    Args
    ----
    sequence    :   Transcript sequence.
    strand      :   '+' (5' to 3') or '-' (3' to 5').
    df_exons    :   DataFrame with the following columns:
                    'exon_id',
                    'exon_number',
                    'start' (reference start position),
                    'end' (reference end position),
                    'size'

    Returns
    -------
    indel_type ('small_insertion' or 'small_deletion')
    variant_sequence,
    variant_start_position (position on reference genome),
    variant_end_position (position on reference genome)
    """
    # Step 1. Reverse sequence if on the minus strand
    if strand == '-':
        sequence = sequence[::-1]

    # Step 2. Fetch the start position of the first exon
    df_first_exon = df_exons.loc[df_exons['exon_number'] == 1,:]
    start = df_first_exon['start'].values[0]

    # Step 3. Randomly select an exon
    exon_number = random.choice(df_exons['exon_number'].values.tolist())
    df_matched = df_exons.loc[df_exons['exon_number'] == exon_number,:]
    exon_size = df_matched['size'].values[0]

    # Step 4. Randomly select a position on the exon
    variant_start_position = random.randint(df_matched['start'].values[0],
                                            df_matched['end'].values[0])
    sequence_position = variant_start_position - start

    # Step 5. Randomize between a small insertion or a deletion
    random_indel = random.choice(['small_insertion', 'small_deletion'])

    # Step 6. Randomly select an INDEL size
    random_indel_size = random.randint(1, min([49, int(exon_size * 0.9)]))

    # Step 7. Apply the INDEL variant
    if random_indel == 'small_insertion':
        variant_sequence = sequence[0:sequence_position] + \
                           randomly_generate_dna_sequence(size=random_indel_size) + \
                           sequence[sequence_position:]
        variant_end_position = variant_start_position
    else:
        variant_sequence = sequence[0:sequence_position] + \
                           sequence[sequence_position + random_indel_size:]
        variant_end_position = variant_start_position + random_indel_size

    if strand == '-':
        variant_sequence = variant_sequence[::-1] # reverse

    # Step 8. Reverse sequence if on the minus strand
    if strand == '-':
        variant_sequence = variant_sequence[::-1]

    return random_indel, \
           variant_sequence, \
           variant_start_position, \
           variant_end_position