import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *
from vstolib.gencode import Gencode


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/leework/Documents/Research/projects/seqdata/references/hg38.fa")

    # Step 2. Load GENCODE
    gencode = Gencode(
        gtf_file="/Users/leework/Documents/Research/projects/seqdata/references/gencode.v41.annotation.gtf",
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    # Step 3. Fetch TP53 sequence and create an alternative 3-prime splice site (chr17:7676201)
    df_transcript_tp53 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000269305']
    df_exons_tp53 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_tp53['transcript_id'].values[0]]
    df_exons_tp53.sort_values(by=['number'], ascending=False, inplace=True) # TP53 is on the reverse strand
    tp53_tumor_sequence = ''
    tp53_normal_sequence = ''
    for _,row in df_exons_tp53.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 4:
            # Create an alternative 3 prime splice site
            reference_position = 7676200
            local_position = len(normal_sequence) - (end - reference_position) - 1
            tumor_sequence = normal_sequence[:local_position+1]
        else:
            tumor_sequence = normal_sequence
        tp53_tumor_sequence = tp53_tumor_sequence + tumor_sequence
        tp53_normal_sequence = tp53_normal_sequence + normal_sequence
    tp53_tumor_sequence = reverse_complement(tp53_tumor_sequence)
    tp53_normal_sequence = reverse_complement(tp53_normal_sequence)

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[tp53_tumor_sequence,tp53_normal_sequence],
        output_fastq_file='../../../test/data/fastq/rna-105-tumor_long-read.fastq',
        num_reads=[1,1],
        stranded=True
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [105],
        'chromosome_1': ['chr17'],
        'position_1': [7676201],
        'strand_1': ['*'],
        'operation_1': ['M'],
        'chromosome_2': ['chr17'],
        'position_2': [7676201],
        'strand_2': ['*'],
        'operation_2': ['M'],
        'variant_size': [''],
        'variant_type': ['A3P'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/rna-105-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/rna-105-tumor_ground_truth.tsv', sep='\t', index=False)
