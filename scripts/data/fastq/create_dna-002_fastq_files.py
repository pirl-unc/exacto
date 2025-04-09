import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/leework/Documents/Research/projects/seqdata/references/hg38.fa")

    # Step 2. Fetch TP53 (chr17:7668421-7687490) sequence
    chromosome = 'chr17'
    start = 7668421
    end = 7687490
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, start - 1, end))

    # Step 3. Create a somatic insertion at position 7674225 (ACGTACGTGGTATGCATGCTGAGACTGAGG)
    reference_position = 7674225
    local_position = length - (end - reference_position) - 1
    sequence_tumor = sequence_normal
    assert sequence_tumor[local_position] == 'C'
    sequence_tumor = sequence_tumor[:local_position+1] + 'ACGTACGTGGTATGCATGCTGAGACTGAGG' + sequence_tumor[local_position+1:]

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[sequence_tumor,sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-002-tumor_long-read.fastq',
        num_reads=[3,3]
    )
    create_fastq_file(
        sequences=[sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-002-normal_long-read.fastq',
        num_reads=[3]
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [2],
        'chromosome_1': ['chr17'],
        'position_1': [7674225],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7674226],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [30],
        'variant_type': ['INS'],
        'variant_sequence': ['ACGTACGTGGTATGCATGCTGAGACTGAGG']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/dna-002-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/dna-002-tumor_ground_truth.tsv', sep='\t', index=False)
