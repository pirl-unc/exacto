import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/leework/Documents/Research/projects/seqdata/references/hg38.fa")

    # Step 2. Fetch ALOX15 (chr17:4630919-4641678) sequence
    alox15_chromosome = 'chr17'
    alox15_start = 4630919
    alox15_end = 4641678
    alox15_length = alox15_end - alox15_start + 1
    alox15_sequence_normal = str(fasta.fetch(alox15_chromosome, alox15_start - 1, alox15_end))

    # Step 3. Fetch TP53 (chr17:7668421-7687490) sequence
    tp53_chromosome = 'chr17'
    tp53_start = 7668421
    tp53_end = 7687490
    tp53_length = tp53_end - tp53_start + 1
    tp53_sequence_normal = str(fasta.fetch(tp53_chromosome, tp53_start - 1, tp53_end))

    # Step 4. Create a somatic translocation (chr17:4637154-TATATACGAGCGTACGTGACTGGTACGTTA-chr17:7674880)
    reference_position_1 = 4637154
    reference_position_2 = 7674880
    local_position_1 = alox15_length - (alox15_end - reference_position_1) - 1
    local_position_2 = tp53_length - (tp53_end - reference_position_2) - 1
    sequence_tumor = alox15_sequence_normal[:local_position_1+1] + 'TATATACGAGCGTACGTGACTGGTACGTTA' + tp53_sequence_normal[local_position_2:]

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[sequence_tumor,alox15_sequence_normal,tp53_sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-005-tumor_long-read.fastq',
        num_reads=[3,3,3]
    )
    create_fastq_file(
        sequences=[alox15_sequence_normal,tp53_sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-005-normal_long-read.fastq',
        num_reads=[3,3]
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [6],
        'chromosome_1': ['chr17'],
        'position_1': [4637154],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7674880],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [''],
        'variant_type': ['BND'],
        'variant_sequence': ['TATATACGAGCGTACGTGACTGGTACGTTA']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/dna-005-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/dna-005-tumor_ground_truth.tsv', sep='\t', index=False)
