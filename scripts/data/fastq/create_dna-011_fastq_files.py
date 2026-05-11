import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/seqdata/references/hg38.fa.gz")

    # Step 2. Fetch TP53 (chr17:7668421-7687490) sequence
    chromosome = 'chr17'
    start = 7668421
    end = 7687490
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create 2 somatic insertions
    # At position 7668421 (AGCGGCGAATATCAGCTACCTCTTAAGATC)
    # At position 7687490 (TATCTCGCGAATTCAGCTACTACTACGGGA)
    sequence_tumor = sequence_normal
    sequence_tumor = 'AGCGGCGAATATCAGCTACCTCTTAAGATC' + sequence_tumor + 'TATCTCGCGAATTCAGCTACTACTACGGGA'

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[sequence_tumor,sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-011-tumor_long-read.fastq',
        num_reads=[15,15]
    )
    create_fastq_file(
        sequences=[sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-011-normal_long-read.fastq',
        num_reads=[15]
    )

    # # Step 5. Create TSV file
    # data = {
    #     'variant_call_id': [15, 16],
    #     'chromosome_1': ['chr17', 'chr17'],
    #     'position_1': [7687490, 7668420],
    #     'strand_1': ['*', '*'],
    #     'operation_1': ['D', 'D'],
    #     'chromosome_2': ['chr17', 'chr17'],
    #     'position_2': [7687491, 7668421],
    #     'strand_2': ['*', '*'],
    #     'operation_2': ['U', 'U'],
    #     'variant_size': [30, 30],
    #     'variant_type': ['INS', 'INS'],
    #     'variant_sequence': ['TATCTCGCGAATTCAGCTACTACTACGGGA', 'AGCGGCGAATATCAGCTACCTCTTAAGATC']
    # }
    # pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/dna-011-tumor_ground_truth.tsv', sep='\t', index=False)
    # pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/dna-011-tumor_ground_truth.tsv', sep='\t', index=False)
