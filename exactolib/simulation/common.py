#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions commonly used
by the simulation module.

Author: Jin Seok (Andy) Lee

Last updated date: July 20, 2022
"""


import random
import pandas as pd
import numpy as np
import logging
from pysam import FastaFile
from typing import Tuple


def randomly_generate_dna_sequence(size: int) -> str:
    """
    Randomly generates a DNA sequence.

    Args
    ----
    size    :   Sequence size (int).

    Returns
    -------
    sequence
    """
    atcg = ['A', 'T', 'C', 'G']
    sequence = [random.choice(atcg) for i in range(0, size)]
    return ''.join(sequence)


def get_exonic_sequence(transcript_id: str,
                        df_gtf: pd.DataFrame,
                        reference_genome_fasta_file: str) -> str:
    """
    Returns the exonic sequence of a transcript.

    Args
    ----
    transcript_id               :   Transcript ID.
    df_gtf                      :   DataFrame with the following columns:
                                    'gene_id',
                                    'gene_name',
                                    'gene_type',
                                    'transcript_id',
                                    'transcript_name',
                                    'transcript_type',
                                    'transcript_support_level',
                                    'exon_id',
                                    'exon_number',
                                    'chrom',
                                    'start',
                                    'end',
                                    'strand',
                                    'level'
    reference_genome_fasta_file :   Reference genome FASTA file.

    Returns
    -------
    exonic sequence
    """
    df_gtf_transcript = df_gtf.loc[df_gtf['transcript_id'] == transcript_id, :]
    df_gtf_transcript = df_gtf_transcript.sort_values(by='exon_number', ascending=True)
    sequence = ''
    fasta_object = FastaFile(reference_genome_fasta_file)
    for index, row in df_gtf_transcript.iterrows():
        curr_chr = row['chrom']
        curr_start = int(row['start'])
        curr_end = int(row['end'])
        exon_sequence = fasta_object.fetch(curr_chr, curr_start, curr_end)
        sequence += exon_sequence
    strand = df_gtf_transcript['strand'].unique()[0]
    if strand == '-':
        sequence = sequence[::-1]
    return sequence


