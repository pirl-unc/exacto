import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz")

    # Step 2. Fetch TP53 (chr17:7668421-7687490) sequence
    chromosome = 'chr17'
    start = 7668421
    end = 7687490
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create a somatic deletion (7673801-7674900)
    reference_position_1 = 7673801
    reference_position_2 = 7674900
    local_position_1 = length - (end - reference_position_1) - 1
    local_position_2 = length - (end - reference_position_2) - 1
    sequence_tumor = sequence_normal
    sequence_tumor = sequence_tumor[:local_position_1] + sequence_tumor[local_position_2+1:]

    # Step 4. Create FASTA files
    create_fasta_file(
        sequences=[sequence_tumor, sequence_normal],
        sequence_names=['scga-mini-dna-005-tumor-1', 'scga-mini-dna-005-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-005-tumor.fasta'
    )
    create_fasta_file(
        sequences=[sequence_normal],
        sequence_names=['scga-mini-dna-005-normal-1'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-005-normal.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [5],
        'chromosome_1': ['chr17'],
        'position_1': [reference_position_1 - 1],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [reference_position_2 + 1],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [reference_position_2 - reference_position_1 + 1],
        'variant_type': ['DEL'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-dna-005-tumor_ground_truth.tsv', sep='\t', index=False)
