import pandas as pd
import pysam
from common import *


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/seqdata/references/hg38.fa.gz")

    # Step 2. Fetch TP53 (chr17:7668421-7687490) sequence
    chromosome = 'chr17'
    start = 7668421
    end = 7687490
    length = end - start + 1
    sequence_normal = str(fasta.fetch(chromosome, 7665000 - 1, 7690000))

    # Step 3. Create a somatic duplication (7674001-7676000)
    sequence_tumor = str(fasta.fetch(chromosome, 7665000 - 1, 7668421)) + \
                     str(fasta.fetch(chromosome, 7668421 - 1, 7687490)) + \
                     str(fasta.fetch(chromosome, 7668421 - 1, 7687490)) + \
                     reverse_complement(str(fasta.fetch(chromosome, 7668421 - 1, 7687490))) + \
                     str(fasta.fetch(chromosome, 7687490, 7690000))

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[sequence_tumor,sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-010-tumor_long-read.fastq',
        num_reads=[15,15]
    )
    create_fastq_file(
        sequences=[sequence_normal],
        output_fastq_file='../../../test/data/fastq/dna-010-normal_long-read.fastq',
        num_reads=[15]
    )

    # # Step 5. Create TSV file
    # data = {
    #     'variant_call_id': [12, 13, 14],
    #     'chromosome_1': ['chr17', 'chr17', 'chr17'],
    #     'position_1': [7687490, 7687490, 7668421],
    #     'strand_1': ['*', '*', '*'],
    #     'operation_1': ['D', 'D', 'U'],
    #     'chromosome_2': ['chr17', 'chr17', 'chr17'],
    #     'position_2': [7668421, 7687490, 7687490],
    #     'strand_2': ['*', '*', '*'],
    #     'operation_2': ['U', 'D', 'U'],
    #     'variant_size': [38140, 0, 38140],
    #     'variant_type': ['INVDUP', 'INVDUP', 'INVDUP'],
    #     'variant_sequence': ['', '', '']
    # }
    # pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/dna-010-tumor_ground_truth.tsv', sep='\t', index=False)
    # pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/dna-010-tumor_ground_truth.tsv', sep='\t', index=False)
