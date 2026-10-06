import os
import sys
import pandas as pd
import pysam
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '')))
from common import *
from vstolib.gencode import Gencode


if __name__ == "__main__":
    # Step 1. Load genome data
    fasta = pysam.FastaFile("/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/hg38_chr17-18.fa.gz")

    # Step 2. Load GENCODE
    gencode = Gencode(
        gtf_file="/Users/ajslee/Documents/Research/projects/project_exacto/exacto/test/data/references/gencode.v41.annotation.chr17-18.gtf.gz",
        version='v41',
        species='human',
        levels=[1,2],
        types=['protein_coding']
    )

    # Step 3. Fetch ASPA
    df_transcript_aspa = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000263080']
    df_exons_aspa = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_aspa['transcript_id'].values[0]]
    df_exons_aspa.sort_values(by=['number'], ascending=True, inplace=True) # ASPA is on the forward strand

    # Step 4. Fetch WSCD1
    df_transcript_wscd1 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000317744']
    df_exons_wscd1 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_wscd1['transcript_id'].values[0]]
    df_exons_wscd1.sort_values(by=['number'], ascending=True, inplace=True) # WSCD1 is on the forward strand

    # Step 5. Create a fusion gene
    def get_exon_sequence(df, exon_number: int):
        df_matched = df[df['number'] == exon_number]
        chromosome = df_matched['chromosome'].values[0]
        start = df_matched['start'].values[0]
        end = df_matched['end'].values[0]
        sequence = str(fasta.fetch(chromosome, start - 1, end))
        return sequence

    aspa_exon_1 = get_exon_sequence(df=df_exons_aspa, exon_number=1)
    aspa_exon_2 = get_exon_sequence(df=df_exons_aspa, exon_number=2)
    aspa_exon_3 = get_exon_sequence(df=df_exons_aspa, exon_number=3)
    aspa_exon_4 = get_exon_sequence(df=df_exons_aspa, exon_number=4)
    aspa_exon_5 = get_exon_sequence(df=df_exons_aspa, exon_number=5)
    aspa_exon_6 = get_exon_sequence(df=df_exons_aspa, exon_number=6)

    wscd1_exon_1 = get_exon_sequence(df=df_exons_wscd1, exon_number=1)
    wscd1_exon_2 = get_exon_sequence(df=df_exons_wscd1, exon_number=2)
    wscd1_exon_3 = get_exon_sequence(df=df_exons_wscd1, exon_number=3)
    wscd1_exon_4 = get_exon_sequence(df=df_exons_wscd1, exon_number=4)
    wscd1_exon_5 = get_exon_sequence(df=df_exons_wscd1, exon_number=5)
    wscd1_exon_6 = get_exon_sequence(df=df_exons_wscd1, exon_number=6)
    wscd1_exon_7 = get_exon_sequence(df=df_exons_wscd1, exon_number=7)
    wscd1_exon_8 = get_exon_sequence(df=df_exons_wscd1, exon_number=8)
    wscd1_exon_9 = get_exon_sequence(df=df_exons_wscd1, exon_number=9)

    tumor_sequence = \
        aspa_exon_1 + \
        aspa_exon_2 + \
        aspa_exon_3 + \
        aspa_exon_4 + \
        wscd1_exon_3 + \
        wscd1_exon_4 + \
        wscd1_exon_5 + \
        wscd1_exon_6 + \
        wscd1_exon_7 + \
        wscd1_exon_8 + \
        wscd1_exon_9
    wscd1_sequence = \
        wscd1_exon_1 + \
        wscd1_exon_2 + \
        wscd1_exon_3 + \
        wscd1_exon_4 + \
        wscd1_exon_5 + \
        wscd1_exon_6 + \
        wscd1_exon_7 + \
        wscd1_exon_8 + \
        wscd1_exon_9
    aspa_sequence = \
        aspa_exon_1 + \
        aspa_exon_2 + \
        aspa_exon_3 + \
        aspa_exon_4 + \
        aspa_exon_5 + \
        aspa_exon_6

    tumor_peptide = longest_orf_peptide(seq=tumor_sequence)
    normal_peptide_1 = longest_orf_peptide(seq=wscd1_sequence)
    normal_peptide_2 = longest_orf_peptide(seq=aspa_sequence)

    # The same transcripts as structures
    def get_exon_structure(df, exon_number: int):
        return TranscriptStructure.exon(fasta, df[df['number'] == exon_number].iloc[0])

    aspa_exon_structures = {number: get_exon_structure(df=df_exons_aspa, exon_number=number) for number in range(1, 7)}
    wscd1_exon_structures = {number: get_exon_structure(df=df_exons_wscd1, exon_number=number) for number in range(1, 10)}
    tumor_structure = aspa_exon_structures[1] + aspa_exon_structures[2] + aspa_exon_structures[3] + aspa_exon_structures[4] + \
        wscd1_exon_structures[3] + wscd1_exon_structures[4] + wscd1_exon_structures[5] + wscd1_exon_structures[6] + wscd1_exon_structures[7] + wscd1_exon_structures[8] + wscd1_exon_structures[9]
    wscd1_structure = wscd1_exon_structures[1] + wscd1_exon_structures[2] + wscd1_exon_structures[3] + wscd1_exon_structures[4] + wscd1_exon_structures[5] + wscd1_exon_structures[6] + wscd1_exon_structures[7] + wscd1_exon_structures[8] + wscd1_exon_structures[9]
    aspa_structure = aspa_exon_structures[1] + aspa_exon_structures[2] + aspa_exon_structures[3] + aspa_exon_structures[4] + aspa_exon_structures[5] + aspa_exon_structures[6]

    # Step 6. Create FASTA files
    create_fasta_file(
        sequences=[tumor_sequence, wscd1_sequence, aspa_sequence],
        sequence_names=['scga-mini-rna-007-tumor-1', 'scga-mini-rna-007-tumor-2', 'scga-mini-rna-007-tumor-3'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-rna-007-tumor.fasta'
    )
    create_fasta_file(
        sequences=[tumor_peptide, normal_peptide_1, normal_peptide_2],
        sequence_names=['scga-mini-pep-007-tumor-1', 'scga-mini-pep-007-tumor-2', 'scga-mini-pep-007-tumor-3'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-pep-007-tumor.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [107],
        'chromosome_1': ['chr17'],
        'position_1': [3489342],
        'strand_1': ['+'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [6087990],
        'strand_2': ['+'],
        'operation_2': ['U'],
        'variant_size': [''],
        'variant_type': ['FUS'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-rna-007-tumor_ground_truth.tsv', sep='\t', index=False)

    # Step 6. Create transcript model alignments TSV file
    write_transcript_model_alignments_ground_truth(
        gencode=gencode,
        transcripts=[
            ('scga-mini-rna-007-tumor-1', tumor_structure, tumor_sequence),
            ('scga-mini-rna-007-tumor-2', wscd1_structure, wscd1_sequence),
            ('scga-mini-rna-007-tumor-3', aspa_structure, aspa_sequence)
        ],
        output_tsv_file='../../../../test/data/simulation/ground_truth/scga-mini-rna-007-tumor_transcript_model_alignments_ground_truth.tsv'
    )
