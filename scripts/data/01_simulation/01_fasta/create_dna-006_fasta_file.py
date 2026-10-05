import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz")

    # Step 2. Fetch WRAP53
    chromosome = 'chr17'
    start = 7687571
    end = 7717000 + 20000
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create a somatic inversion (7701200-7717000)
    sequence_tumor = str(fasta.fetch(chromosome, start - 1, 7701200)) + \
                     reverse_complement(str(fasta.fetch(chromosome, 7701200, 7717000))) + \
                     str(fasta.fetch(chromosome, 7717000, 7717000 + 20000))

    # Step 4. Create FASTA files
    create_fasta_file(
        sequences=[sequence_tumor, sequence_normal],
        sequence_names=['scga-mini-dna-006-tumor-1', 'scga-mini-dna-006-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-006-tumor.fasta'
    )
    create_fasta_file(
        sequences=[sequence_normal],
        sequence_names=['scga-mini-dna-006-normal-1'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-dna-006-normal.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [6,7],
        'chromosome_1': ['chr17','chr17'],
        'position_1': [7701200,7701200+1],
        'strand_1': ['*','*'],
        'operation_1': ['D','U'],
        'chromosome_2': ['chr17','chr17'],
        'position_2': [7717000,7717000+1],
        'strand_2': ['*','*'],
        'operation_2': ['D','U'],
        'variant_size': ['',''],
        'variant_type': ['BND','BND'],
        'variant_sequence': ['','']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-dna-006-tumor_ground_truth.tsv', sep='\t', index=False)
