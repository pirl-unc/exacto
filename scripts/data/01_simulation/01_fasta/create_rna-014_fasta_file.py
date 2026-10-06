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

    # Step 3. Fetch TP53 sequence and create a somatic insertion in the middle of exon 11 of TP53 sequence (TTAACTGAGTCTCAAAAAAATAAA).
    # The insertion is written on the chr17 + strand and TP53 is on the reverse strand, so the mRNA carries the reverse complement
    # (TTTATTTTTTTGAGACTCAGTTAA), which reads FIFLRLS and then TAA in the TP53 frame - a premature stop after residue 380.
    df_transcript_tp53 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000269305']
    df_exons_tp53 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_tp53['transcript_id'].values[0]]
    df_exons_tp53.sort_values(by=['number'], ascending=False, inplace=True) # TP53 is on the reverse strand

    def get_exon_sequence(exon_number: int):
        df_matched = df_exons_tp53[df_exons_tp53['number'] == exon_number]
        chromosome = df_matched['chromosome'].values[0]
        start = df_matched['start'].values[0]
        end = df_matched['end'].values[0]
        sequence = str(fasta.fetch(chromosome, start - 1, end))
        return sequence

    exon_1 = get_exon_sequence(exon_number=1)
    exon_2 = get_exon_sequence(exon_number=2)
    exon_3 = get_exon_sequence(exon_number=3)
    exon_4 = get_exon_sequence(exon_number=4)
    exon_5 = get_exon_sequence(exon_number=5)
    exon_6 = get_exon_sequence(exon_number=6)
    exon_7 = get_exon_sequence(exon_number=7)
    exon_8 = get_exon_sequence(exon_number=8)
    exon_9 = get_exon_sequence(exon_number=9)
    exon_10 = get_exon_sequence(exon_number=10)
    exon_11 = get_exon_sequence(exon_number=11)

    # Tumor
    start = df_exons_tp53.loc[df_exons_tp53['number'] == 11, 'start'].values[0]
    end = df_exons_tp53.loc[df_exons_tp53['number'] == 11, 'end'].values[0]
    exon_11_tumor = \
        str(fasta.fetch('chr17', start - 1, 7669650)) + \
        "TTAACTGAGTCTCAAAAAAATAAA" + \
        str(fasta.fetch('chr17', 7669651 - 1, end))
    tp53_tumor_sequence = \
        exon_11_tumor + \
        exon_10 + \
        exon_9 + \
        exon_8 + \
        exon_7 + \
        exon_6 + \
        exon_5 + \
        exon_4 + \
        exon_3 + \
        exon_2 + \
        exon_1
    tp53_tumor_sequence = reverse_complement(tp53_tumor_sequence)

    # Normal
    tp53_normal_sequence = \
        exon_11 + \
        exon_10 + \
        exon_9 + \
        exon_8 + \
        exon_7 + \
        exon_6 + \
        exon_5 + \
        exon_4 + \
        exon_3 + \
        exon_2 + \
        exon_1
    tp53_normal_sequence = reverse_complement(tp53_normal_sequence)

    tp53_tumor_peptide = longest_orf_peptide(seq=tp53_tumor_sequence)
    tp53_normal_peptide = longest_orf_peptide(seq=tp53_normal_sequence)

    # The same transcripts as structures
    def get_exon_structure(exon_number: int):
        return TranscriptStructure.exon(fasta, df_exons_tp53[df_exons_tp53['number'] == exon_number].iloc[0])

    exon_11_row = df_exons_tp53[df_exons_tp53['number'] == 11].iloc[0]
    exon_11_tumor_structure = \
        TranscriptStructure.exon(fasta, exon_11_row, end=7669650) + \
        TranscriptStructure.insertion("TTAACTGAGTCTCAAAAAAATAAA") + \
        TranscriptStructure.exon(fasta, exon_11_row, start=7669651)
    tp53_tumor_structure = (exon_11_tumor_structure + get_exon_structure(10) + get_exon_structure(9) + get_exon_structure(8) + get_exon_structure(7) + get_exon_structure(6) + get_exon_structure(5) + get_exon_structure(4) + get_exon_structure(3) + get_exon_structure(2) + get_exon_structure(1)).reverse_complement()
    tp53_normal_structure = (get_exon_structure(11) + get_exon_structure(10) + get_exon_structure(9) + get_exon_structure(8) + get_exon_structure(7) + get_exon_structure(6) + get_exon_structure(5) + get_exon_structure(4) + get_exon_structure(3) + get_exon_structure(2) + get_exon_structure(1)).reverse_complement()

    # Step 4. Create FASTA files
    create_fasta_file(
        sequences=[tp53_tumor_sequence, tp53_normal_sequence],
        sequence_names=['scga-mini-rna-014-tumor-1', 'scga-mini-rna-014-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-rna-014-tumor.fasta'
    )
    create_fasta_file(
        sequences=[tp53_tumor_peptide, tp53_normal_peptide],
        sequence_names=['scga-mini-pep-014-tumor-1', 'scga-mini-pep-014-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-pep-014-tumor.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [115],
        'chromosome_1': ['chr17'],
        'position_1': [7669650],
        'strand_1': ['-'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7669651],
        'strand_2': ['-'],
        'operation_2': ['U'],
        'variant_size': [24],
        'variant_type': ['INS'],
        'variant_sequence': ['TTTATTTTTTTGAGACTCAGTTAA']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-rna-014-tumor_ground_truth.tsv', sep='\t', index=False)

    # Step 6. Create transcript model alignments TSV file
    write_transcript_model_alignments_ground_truth(
        gencode=gencode,
        transcripts=[
            ('scga-mini-rna-014-tumor-1', tp53_tumor_structure, tp53_tumor_sequence),
            ('scga-mini-rna-014-tumor-2', tp53_normal_structure, tp53_normal_sequence)
        ],
        output_tsv_file='../../../../test/data/simulation/ground_truth/scga-mini-rna-014-tumor_transcript_model_alignments_ground_truth.tsv'
    )
