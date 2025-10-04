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

    # Step 3. Fetch METTL16 sequence
    df_transcript_mettl16 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000263092']
    df_exons_mettl16 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_mettl16['transcript_id'].values[0]]
    df_exons_mettl16.sort_values(by=['number'], ascending=False, inplace=True) # METTL16 is on the reverse strand
    mettl16_tumor_sequence = ''
    mettl16_normal_sequence = ''
    for _,row in df_exons_mettl16.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 6:
            # Create a fusion (chr17:2464300-4433940)
            reference_position = 2464300
            local_position = len(normal_sequence) - (end - reference_position) - 1
            mettl16_tumor_sequence = mettl16_tumor_sequence + normal_sequence[:local_position+1]
        if exon_number < 6:
            mettl16_tumor_sequence = mettl16_tumor_sequence + normal_sequence
        mettl16_normal_sequence = mettl16_normal_sequence + normal_sequence
    mettl16_tumor_sequence = reverse_complement(mettl16_tumor_sequence)
    mettl16_normal_sequence = reverse_complement(mettl16_normal_sequence)

    # Step 4. Fetch SPNS3 sequence
    df_transcript_spns3 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000355530']
    df_exons_spns3 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_spns3['transcript_id'].values[0]]
    df_exons_spns3.sort_values(by=['number'], ascending=True, inplace=True) # ITGAE is on the forward strand
    spns3_tumor_sequence = ''
    spns3_normal_sequence = ''
    for _,row in df_exons_spns3.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 8:
            # Create a fusion (chr17:4453100-7603800)
            reference_position = 4453100
            local_position = len(normal_sequence) - (end - reference_position) - 1
            spns3_tumor_sequence = spns3_tumor_sequence + normal_sequence[:local_position+1]
        if exon_number < 8:
            spns3_tumor_sequence = spns3_tumor_sequence + normal_sequence
        spns3_normal_sequence = spns3_normal_sequence + normal_sequence

    # Step 5. Fetch FXR2 sequence
    df_transcript_fxr2 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000250113']
    df_exons_fxr2 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_fxr2['transcript_id'].values[0]]
    df_exons_fxr2.sort_values(by=['number'], ascending=False, inplace=True) # FXR2 is on the reverse strand
    fxr2_tumor_sequence = ''
    fxr2_normal_sequence = ''
    for _,row in df_exons_fxr2.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 5:
            # Create a fusion (chr17:7603800-7591230)
            reference_position = 7603800
            local_position = len(normal_sequence) - (end - reference_position) - 1
            fxr2_tumor_sequence = fxr2_tumor_sequence + normal_sequence[:local_position]
        if exon_number > 5:
            fxr2_tumor_sequence = fxr2_tumor_sequence + normal_sequence
        fxr2_normal_sequence = fxr2_normal_sequence + normal_sequence
    fxr2_tumor_sequence = reverse_complement(fxr2_tumor_sequence)
    fxr2_normal_sequence = reverse_complement(fxr2_normal_sequence)

    # Step 6. Create FASTQ files
    tumor_sequence = mettl16_tumor_sequence + spns3_tumor_sequence + fxr2_tumor_sequence
    create_fastq_file(
        sequences=[tumor_sequence, mettl16_normal_sequence, spns3_normal_sequence, fxr2_normal_sequence],
        output_fastq_file='../../../test/data/fastq/rna-110-tumor_long-read.fastq',
        num_reads=[1,1,1,1],
        stranded=True
    )

    # Step 7. Create TSV file
    data = {
        'variant_call_id': [1100,1101],
        'chromosome_1': ['chr17','chr17'],
        'position_1': [2464300, 4453100],
        'strand_1': ['*', '*'],
        'operation_1': ['U', 'D'],
        'chromosome_2': ['chr17', 'chr17'],
        'position_2': [4433940, 7603800],
        'strand_2': ['*', '*'],
        'operation_2': ['U', 'D'],
        'variant_size': ['', ''],
        'variant_type': ['FUS', 'FUS'],
        'variant_sequence': ['', '']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/rna-110-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/rna-110-tumor_ground_truth.tsv', sep='\t', index=False)

