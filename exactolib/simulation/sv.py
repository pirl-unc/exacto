#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
simulating transcript-level structural variants.

Author: Jin Seok (Andy) Lee

Last updated date: July 20, 2022
"""


import random
import pandas as pd
import numpy as np
import logging
from typing import Tuple
from exactolib.simulation.common import *


def randomly_generate_deletion(transcript_id: str,
                               sequence: str,
                               strand: str,
                               df_gtf_transcript: pd.DataFrame,
                               df_trasncript_sequences: pd.DataFrame) -> Tuple[str, str, int, int, int, int]:
    """
    Randomly generates a deletion (>= 50-bp).

    Args
    ----
    transcript_id               :   Transcript ID (Ensemble transcript ID).
    sequence                    :   Transcript sequence.
    strand                      :   '+' (5' to 3') or '-' (3' to 5').
    df_gtf_transcript           :   This DataFrame should only contain all the exons
                                    of one transcript. DataFrame with the following columns:
                                    'transcript_id',
                                    'exon_id',
                                    'exon_number',
                                    'chrom',
                                    'start' (reference start position),
                                    'end' (reference end position),
                                    'size'.
    df_trasncript_sequences     :   DataFrame with the following columns:
                                    'ensembl_transcript_id',
                                    'ensembl_gene_id',
                                    'utr5_start',
                                    'utr5_end',
                                    'cds_start',
                                    'cds_end',
                                    'utr3_start',
                                    'utr3_end',
                                    'sequence'

    Returns
    -------
    deletion_type ('deletion_type_1',
                   'deletion_type_2',
                   'deletion_type_3',
                   'deletion_type_4', or
                   'deletion_type_5')
    variant_sequence,
    variant_start_position (position on reference genome),
    variant_end_position (position on reference genome),
    deletion_size,
    exon_number
    """
    # Step 1. Reverse sequence if on the minus strand
    if strand == '-':
        sequence = sequence[::-1]

    # Step 2. Fetch the start position of the first exon
    df_first_exon = df_gtf_transcript.loc[df_gtf_transcript['exon_number'] == 1,:]
    start = df_first_exon['start'].values[0]

    # Step 3. Randomly select an exon
    df_exons_large = df_gtf_transcript.loc[df_gtf_transcript['size'] > 50, :]
    exon_number = random.choice(df_exons_large['exon_number'].values.tolist())
    df_matched = df_exons_large.loc[df_exons_large['exon_number'] == exon_number, :]
    exon_size = df_matched['size'].values[0]

    # Step 4. Randomly select a deletion type
    if (1 == exon_number) or (len(df_gtf_transcript) == exon_number):
        # If either the first or the last exon was selected,
        # then always simulate deletion type 4
        deletion_type = 'deletion_type_4'
    else:
        deletion_types = ['deletion_type_1',
                          'deletion_type_2',
                          'deletion_type_3',
                          'deletion_type_5']
        deletion_type = random.choice(deletion_types)

    # Step 5. Apply the deletion type
    if deletion_type == 'deletion_type_1': # deletion in the middle of an exon
        # Randomly select deletion start position
        variant_start_position = random.randint(df_matched['start'].values[0] + 1,
                                                df_matched['end'].values[0] - 1)
        sequence_start_position = variant_start_position - start

        # Randomly select a deletion size
        deletion_size = random.randint(50, int(exon_size * 0.9))
        variant_end_position = variant_start_position + deletion_size
        sequence_end_position = sequence_start_position + deletion_size

        # Apply the deletion
        variant_sequence = sequence[0:sequence_start_position] + \
                            sequence[sequence_end_position:]
    elif deletion_type == 'deletion_type_2': # alternative 5' splice site
        variant_start_position = df_matched['start'].values[0]
        sequence_start_position = 0

        # Randomly select a deletion size
        deletion_size = random.randint(50, int(exon_size * 0.9))
        variant_end_position = variant_start_position + deletion_size
        sequence_end_position = sequence_start_position + deletion_size

        # Apply the deletion
        variant_sequence = sequence[0:sequence_start_position] + \
                            sequence[sequence_end_position:]
    elif deletion_type == 'deletion_type_3': # alternative 3' splice site
        # Randomly select a deletion size
        deletion_size = random.randint(50, int(exon_size * 0.9))
        variant_end_position = df_matched['end'].values[0]
        variant_start_position = variant_end_position - deletion_size
        sequence_start_position = variant_start_position - start
        sequence_end_position = sequence_start_position + deletion_size

        # Apply the deletion
        variant_sequence = sequence[0:sequence_start_position] + \
                            sequence[sequence_end_position:]
    elif deletion_type == 'deletion_type_4': # 5' of first exon or 3' of last exon
        if exon_number == 1: # 5' deletion
            # Randomly select a deletion size
            deletion_size = random.randint(50, int(exon_size * 0.9))
            variant_end_position = df_matched['start'].values[0] + deletion_size
            variant_start_position = df_matched['start'].values[0]
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position + deletion_size

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                                sequence[sequence_end_position:]
        else: # 3' deletion
            # Randomly select a deletion size
            deletion_size = random.randint(50, int(exon_size * 0.9))
            variant_end_position = df_matched['end'].values[0]
            variant_start_position = df_matched['end'].values[0] - deletion_size
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position + deletion_size

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               sequence[sequence_end_position:]
    elif deletion_type == 'deletion_type_5': # exon skipping
        # Ensure that the exon skipping does not lead to`another known isoform
        curr_transcript_sequences = df_trasncript_sequences.loc[df_trasncript_sequences['ensembl_transcript_id'] == transcript_id, 'sequence'].values.tolist()
        while True:
            # Generate sequence without the currently selected exon
            df_gtf_transcript_exon_skipped = df_gtf_transcript.loc[df_gtf_transcript['exon_number'] != exon_number,:]
            df_gtf_transcript_exon_skipped.sort_values(by='exon_number', ascending=True, inplace=True)
            curr_sequence = ''
            for index, row in df_gtf_transcript_exon_skipped.iterrows():
                curr_sequence += row['']

        # Exon skipping
        deletion_size = df_matched['end'].values[0] - df_matched['start'].values[0] + 1
        variant_end_position = df_matched['end'].values[0] + 1
        variant_start_position = df_matched['start'].values[0] - 1
        sequence_start_position = variant_start_position - start
        sequence_end_position = sequence_start_position + deletion_size

        # Apply the deletion
        variant_sequence = sequence[0:sequence_start_position] + \
                           sequence[sequence_end_position:]
    else:
        print("unknown deletion type selected")

    # Step 6. Reverse sequence if on the minus strand
    if strand == '-':
        variant_sequence = variant_sequence[::-1]

    return deletion_type, \
           variant_sequence, \
           variant_start_position,\
           variant_end_position, \
           deletion_size, \
           exon_number


def randomly_generate_insertion(sequence: str,
                                strand: str,
                                df_exons: pd.DataFrame,
                                max_size=500) -> Tuple[str, str, ]:
    """
    Randomly generates an insertion (>= 50-bp).

    Args
    ----
    sequence    :   nucleotide sequence.
    strand      :   '+' (5' to 3') or '-' (3' to 5').
    df_exons    :   DataFrame with the following columns
                    'exon_id',
                    'exon_number',
                    'start' (reference start position),
                    'end' (reference end position),
                    'size'
    max_size    :   maximum insertion size (default: 1,000).

    Returns
    -------
    variant_sequence,
    variant_start_position (position on reference genome),
    variant_end_position (position on reference genome),
    insertion_sequence,
    insertion_size,
    exon_number
    """
    # Step 1. Reverse sequence if on the minus strand
    if strand == '-':
        sequence = sequence[::-1]

    # Step 2. Fetch the start position of the first exon
    df_first_exon = df_exons.loc[df_exons['exon_number'] == 1,:]
    start = df_first_exon['start'].values[0]

    # Step 3. Randomly select an exon
    exon_number = random.choice(df_exons['exon_number'].values.tolist())
    df_matched = df_exons.loc[df_exons['exon_number'] == exon_number, :]
    exon_size = df_matched['size'].values[0]

    # Step 4. Randomly select an insertion type
    if (1 == exon_number) or (len(df_exons) == exon_number):
        # If either the first or the last exon was selected, then always
        # simulate insertion type 3
        random_insertion_type = 'insertion_type_3'
    else:
        insertion_types = ['insertion_type_1',
                           'insertion_type_2',
                           'insertion_type_4']
        random_insertion_type = random.choice(insertion_types)

    # Step 5. Apply the insertion type
    if random_insertion_type == 'insertion_type_1':
        # Randomly select insertion start position
        variant_start_position = random.randint(df_matched['start'].values[0] + 1,
                                                df_matched['end'].values[0] - 1)
        variant_end_position = variant_start_position
        sequence_start_position = variant_start_position - start
        sequence_end_position = sequence_start_position

        # Randomly select an insertion size
        insertion_size = random.randint(50, max_size)
        insertion_sequence = randomly_generate_dna_sequence(size=insertion_size)

        # Apply the deletion
        variant_sequence = sequence[0:sequence_start_position] + \
                            insertion_sequence + \
                            sequence[sequence_end_position:]
    elif random_insertion_type == 'insertion_type_2':
        # Randomly select an insertion size
        insertion_size = random.randint(50, max_size)
        insertion_sequence = randomly_generate_dna_sequence(size=insertion_size)

        # Randomize between 5' and 3' insertion
        insertion_site = random.choice(['5', '3'])
        if insertion_site == '5':
            variant_start_position = df_matched['start'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]
        else: # 3'
            variant_start_position = df_matched['end'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]
    elif random_insertion_type == 'insertion_type_3':
        # Randomly select an insertion size
        insertion_size = random.randint(50, max_size)
        insertion_sequence = randomly_generate_dna_sequence(size=insertion_size)

        # Randomize between 5' and 3' insertion
        if exon_number == 1:
            variant_start_position = df_matched['start'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]
        else: # last exon
            variant_start_position = df_matched['end'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]
    elif random_insertion_type == 'insertion_type_4': # cryptic exon
        # Randomly select an insertion size
        insertion_size = random.randint(50, max_size)
        insertion_sequence = randomly_generate_dna_sequence(size=insertion_size)

        # Insert to the 3' end
        variant_start_position = df_matched['end'].values[0]
        variant_end_position = variant_start_position
        sequence_start_position = variant_start_position - start
        sequence_end_position = sequence_start_position

        # Randomize between 5' and 3' insertion
        if exon_number == 1:
            variant_start_position = df_matched['start'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]
        else: # last exon
            variant_start_position = df_matched['end'].values[0]
            variant_end_position = variant_start_position
            sequence_start_position = variant_start_position - start
            sequence_end_position = sequence_start_position

            # Apply the deletion
            variant_sequence = sequence[0:sequence_start_position] + \
                               insertion_sequence + \
                               sequence[sequence_end_position:]

    # Step 6. Reverse sequence if on the minus strand
    if strand == '-':
        variant_sequence = variant_sequence[::-1]

    return variant_sequence, variant_start_position, variant_end_position, insertion_sequence, insertion_size, exon_number


def randomly_generate_duplication(sequence, strand, df_exons, max_size=500):
    """
    Randomly generates an insertion (>= 50-bp).

    Args
    ----
    sequence    :   nucleotide sequence.
    strand      :   '+' (5' to 3') or '-' (3' to 5').
    df_exons    :   DataFrame with the following columns
                    'exon_id',
                    'exon_number',
                    'start' (reference start position),
                    'end' (reference end position),
                    'size'
    max_size    :   maximum insertion size (default: 1,000).

    Returns
    -------
    variant_sequence,
    variant_start_position (position on reference genome),
    variant_end_position (position on reference genome),
    insertion_sequence,
    insertion_size,
    exon_number
    """
    # Step 1. Reverse sequence if on the minus strand
    if strand == '-':
        sequence = sequence[::-1]

    # Step 2. Fetch the start position of the first exon
    df_first_exon = df_exons.loc[df_exons['exon_number'] == 1,:]
    start = df_first_exon['start'].values[0]
