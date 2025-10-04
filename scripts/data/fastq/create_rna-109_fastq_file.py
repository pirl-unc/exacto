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

    # Step 3. Fetch TRARG1 sequence
    df_transcript_trarg1 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000333813']
    df_exons_trarg1 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_trarg1['transcript_id'].values[0]]
    df_exons_trarg1.sort_values(by=['number'], ascending=True, inplace=True) # TRARG1 is on the forward strand
    trarg1_tumor_sequence = ''
    trarg1_normal_sequence = ''
    for _,row in df_exons_trarg1.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 2:
            # Create a fusion (chr17:1295600-3801188)
            reference_position = 1295600
            local_position = len(normal_sequence) - (end - reference_position) - 1
            trarg1_tumor_sequence = trarg1_tumor_sequence + normal_sequence[:local_position+1]
        if exon_number < 2:
            trarg1_tumor_sequence = trarg1_tumor_sequence + normal_sequence
        trarg1_normal_sequence = trarg1_normal_sequence + normal_sequence

    # Step 4. Fetch ITGAE sequence
    df_transcript_itgae = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000263087']
    df_exons_itgae = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_itgae['transcript_id'].values[0]]
    df_exons_itgae.sort_values(by=['number'], ascending=False, inplace=True) # ITGAE is on the reverse strand
    itgae_tumor_sequence = ''
    itgae_normal_sequence = ''
    for _,row in df_exons_itgae.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 6:
            # Create a fusion (chr17:3761100-7727200)
            reference_position = 3761100
            local_position = len(normal_sequence) - (end - reference_position) - 1
            assert normal_sequence[local_position] == 'C'
            itgae_tumor_sequence = normal_sequence[local_position+1:]
        if exon_number < 6:
            itgae_tumor_sequence = itgae_tumor_sequence + normal_sequence
        itgae_normal_sequence = itgae_normal_sequence + normal_sequence
    itgae_tumor_sequence = reverse_complement(itgae_tumor_sequence)
    itgae_normal_sequence = reverse_complement(itgae_normal_sequence)

    # Step 5. Fetch DNAH2 sequence
    df_transcript_dnah2 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000570791']
    df_exons_dnah2 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_dnah2['transcript_id'].values[0]]
    df_exons_dnah2.sort_values(by=['number'], ascending=True, inplace=True) # DNAH2 is on the forward strand
    dnah2_tumor_sequence = ''
    dnah2_normal_sequence = ''
    for _,row in df_exons_dnah2.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 4:
            # Create a fusion (chr17:7727200-7743901)
            reference_position = 7727200
            local_position = len(normal_sequence) - (end - reference_position) - 1
            dnah2_tumor_sequence = normal_sequence[local_position:]
        if exon_number > 4:
            dnah2_tumor_sequence = dnah2_tumor_sequence + normal_sequence
        dnah2_normal_sequence = dnah2_normal_sequence + normal_sequence

    # Step 6. Create FASTQ files
    tumor_sequence = trarg1_tumor_sequence + itgae_tumor_sequence + dnah2_tumor_sequence
    create_fastq_file(
        sequences=[tumor_sequence, trarg1_normal_sequence, itgae_normal_sequence, dnah2_normal_sequence],
        output_fastq_file='../../../test/data/fastq/rna-109-tumor_long-read.fastq',
        num_reads=[1,1,1,1],
        stranded=True
    )

    # Step 7. Create TSV file
    data = {
        'variant_call_id': [1090,1091],
        'chromosome_1': ['chr17','chr17'],
        'position_1': [1295600, 3761100],
        'strand_1': ['*', '*'],
        'operation_1': ['D', 'U'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [3801188, 7727200],
        'strand_2': ['*', '*'],
        'operation_2': ['D', 'U'],
        'variant_size': ['', ''],
        'variant_type': ['FUS', 'FUS'],
        'variant_sequence': ['', '']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/rna-109-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/rna-109-tumor_ground_truth.tsv', sep='\t', index=False)

