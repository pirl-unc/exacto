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

    # Step 3. Fetch RPA1 sequence
    df_transcript_rpa1 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000254719']
    df_exons_rpa1 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_rpa1['transcript_id'].values[0]]
    df_exons_rpa1.sort_values(by=['number'], ascending=True, inplace=True) # RPA1 is on the forward strand
    rpa1_tumor_sequence = ''
    rpa1_normal_sequence = ''
    for _,row in df_exons_rpa1.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number >= 2 and exon_number <= 4:
            rpa1_tumor_sequence = rpa1_tumor_sequence + normal_sequence
        rpa1_normal_sequence = rpa1_normal_sequence + normal_sequence

    # Step 4. Create FASTQ files
    tumor_sequence = rpa1_tumor_sequence + rpa1_tumor_sequence
    create_fastq_file(
        sequences=[tumor_sequence, rpa1_normal_sequence],
        output_fastq_file='../../../test/data/fastq/rna-111-tumor_long-read.fastq',
        num_reads=[1,1],
        stranded=True
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [111],
        'chromosome_1': ['chr17'],
        'position_1': [1844686],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [1842803],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [''],
        'variant_type': ['CIR'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/rna-111-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/rna-111-tumor_ground_truth.tsv', sep='\t', index=False)

