import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz")

    # Step 2. Fetch ASPA (chr17:3475997-3503405) sequence
    chromosome = 'chr17'
    start = 3475997
    end = 3503405
    length = end - start + 1
    aspa_sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Fetch WSCD1 (chr17:6069106-6124427) sequence
    chromosome = 'chr17'
    start = 6069106
    end = 6124427
    length = end - start + 1
    wscd1_sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 4. Create a translocation / deletion (chr17:3491600-6085000)
    sequence_tumor = str(fasta.fetch(chromosome, 3475997 - 1, 3491600)) + str(fasta.fetch(chromosome, 6085000, 6124427))

    # Step 5. Create FASTA files
    create_fasta_file(
        sequences=[aspa_sequence_normal, wscd1_sequence_normal, sequence_tumor],
        sequence_names=['scga-mini-dna-007-tumor-1', 'scga-mini-dna-007-tumor-2', 'scga-mini-dna-007-tumor-3'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-007-tumor.fasta'
    )
    create_fasta_file(
        sequences=[aspa_sequence_normal, wscd1_sequence_normal],
        sequence_names=['scga-mini-dna-007-normal-1', 'scga-mini-dna-007-normal-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-007-normal.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [8],
        'chromosome_1': ['chr17'],
        'position_1': [3491600],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [6085000+1],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [''],
        'variant_type': ['BND'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-dna-007-tumor_ground_truth.tsv', sep='\t', index=False)
