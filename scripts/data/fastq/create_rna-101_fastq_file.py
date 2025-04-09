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

    # Step 3. Fetch TP53 sequence and create a somatic insertion at position 7676400 (GGGGGTTTTT)
    df_transcript_tp53 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000269305']
    df_exons_tp53 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_tp53['transcript_id'].values[0]]
    df_exons_tp53.sort_values(by=['number'], ascending=False, inplace=True) # TP53 is on the reverse strand
    tp53_tumor_sequence = ''
    tp53_normal_sequence = ''
    for _,row in df_exons_tp53.iterrows():
        print(row)
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        if exon_number == 3:
            # Create a SNV at position 7676400
            reference_position = 7676400
            local_position = len(normal_sequence) - (end - reference_position) - 1
            assert normal_sequence[local_position] == 'A'
            tumor_sequence = normal_sequence[:local_position+1] + 'GGGGGTTTTT' + normal_sequence[local_position+1:]
        else:
            tumor_sequence = normal_sequence
        tp53_tumor_sequence = tp53_tumor_sequence + tumor_sequence
        tp53_normal_sequence = tp53_normal_sequence + normal_sequence
    tp53_tumor_sequence = reverse_complement(tp53_tumor_sequence)
    tp53_normal_sequence = reverse_complement(tp53_normal_sequence)

    # Step 4. Create FASTQ files
    create_fastq_file(
        sequences=[tp53_tumor_sequence,tp53_normal_sequence],
        output_fastq_file='../../../test/data/fastq/rna-101-tumor_long-read.fastq',
        num_reads=[1,1],
        stranded=True
    )

    # Step 5. Create TSV file
    data = {
        'variant_call_id': [101],
        'chromosome_1': ['chr17'],
        'position_1': [7676400],
        'strand_1': ['*'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7676401],
        'strand_2': ['*'],
        'operation_2': ['U'],
        'variant_size': [10],
        'variant_type': ['INS'],
        'variant_sequence': ['GGGGGTTTTT']
    }
    pd.DataFrame(data).to_csv('../../../test/data/tsv/ground_truth/rna-101-tumor_ground_truth.tsv', sep='\t', index=False)
    pd.DataFrame(data).to_csv('/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/tsv/ground_truth/rna-101-tumor_ground_truth.tsv', sep='\t', index=False)
