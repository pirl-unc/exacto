import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/leework/Documents/Research/projects/seqdata/references/hg38.fa")

    # Step 2. Fetch AKAIN1 (chr18:5142911-5197257) sequence
    akain1_chromosome = 'chr18'
    akain1_start = 5142911
    akain1_end = 5197257
    akain1_length = akain1_end - akain1_start + 1
    akain1_sequence_normal = str(fasta.fetch(akain1_chromosome, akain1_start - 1, akain1_end))

    # Step 3. Fetch TP53 (chr17:7668421-7687490) sequence
    tp53_chromosome = 'chr17'
    tp53_start = 7668421
    tp53_end = 7687490
    tp53_length = tp53_end - tp53_start + 1
    tp53_sequence_normal = str(fasta.fetch(tp53_chromosome, tp53_start - 1, tp53_end))

    # Step 4. Create a somatic translocation (chr17:7676155-chr18:5170100)
    reference_position_1 = 7676155
    reference_position_2 = 5170100
    local_position_1 = tp53_length - (tp53_end - reference_position_1) - 1
    local_position_2 = akain1_length - (akain1_end - reference_position_2) - 1
    sequence_tumor = tp53_sequence_normal[:local_position_1+1] + akain1_sequence_normal[local_position_2:]

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[sequence_tumor,akain1_sequence_normal,tp53_sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-006-tumor_long-read.fastq',
        num_reads=[3,3,3]
    )
    create_fastq_file(
        sequences=[akain1_sequence_normal,tp53_sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-006-normal_long-read.fastq',
        num_reads=[3,3]
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [7],
        'chromosome_1': ['chr17'],
        'position_1': [7676155],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr18'],
        'position_2': [5170100],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [''],
        'variant_type': ['BND'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/dna-006-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/dna-006-tumor_ground_truth.tsv', sep='\t', index=False)
