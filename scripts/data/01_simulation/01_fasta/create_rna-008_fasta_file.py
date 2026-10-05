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

    # Step 3. Fetch TP53 sequence and create a cryptic exon (7672205:7672484)
    df_transcript_tp53 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000269305']
    df_exons_tp53 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_tp53['transcript_id'].values[0]]
    df_exons_tp53.sort_values(by=['number'], ascending=False, inplace=True) # TP53 is on the reverse strand
    tp53_tumor_sequence = ''
    tp53_normal_sequence = ''
    tp53_tumor_structure = TranscriptStructure()
    tp53_normal_structure = TranscriptStructure()
    for _,row in df_exons_tp53.iterrows():
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        exon_number = row['number']
        normal_sequence = str(fasta.fetch(chromosome, start - 1, end))
        normal_structure = TranscriptStructure.exon(fasta, row)
        if exon_number == 9:
            # Create a cryptic exon (7672205:7672484)
            cryptic_exon = str(fasta.fetch(chromosome, 7672205 - 1, 7672484))
            tumor_sequence = cryptic_exon + normal_sequence
            tumor_structure = TranscriptStructure.segment(fasta, chromosome, 7672205, 7672484, row['gene_id'], row['transcript_id']) + \
                              normal_structure
        else:
            tumor_sequence = normal_sequence
            tumor_structure = normal_structure
        tp53_tumor_sequence = tp53_tumor_sequence + tumor_sequence
        tp53_normal_sequence = tp53_normal_sequence + normal_sequence
        tp53_tumor_structure = tp53_tumor_structure + tumor_structure
        tp53_normal_structure = tp53_normal_structure + normal_structure
    tp53_tumor_sequence = reverse_complement(tp53_tumor_sequence)
    tp53_normal_sequence = reverse_complement(tp53_normal_sequence)
    tp53_tumor_structure = tp53_tumor_structure.reverse_complement()
    tp53_normal_structure = tp53_normal_structure.reverse_complement()

    tp53_tumor_peptide = longest_orf_peptide(seq=tp53_tumor_sequence)
    tp53_normal_peptide = longest_orf_peptide(seq=tp53_normal_sequence)

    # Step 4. Create FASTA files
    create_fasta_file(
        sequences=[tp53_tumor_sequence, tp53_normal_sequence],
        sequence_names=['scga-mini-rna-008-tumor-1', 'scga-mini-rna-008-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-rna-008-tumor.fasta'
    )
    create_fasta_file(
        sequences=[tp53_tumor_peptide, tp53_normal_peptide],
        sequence_names=['scga-mini-pep-008-tumor-1', 'scga-mini-pep-008-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-pep-008-tumor.fasta'
    )

    # Step 5. Create TSV file
    data = {
        'variant_id': [108],
        'chromosome_1': ['chr17'],
        'position_1': [7672205],
        'strand_1': ['-'],
        'operation_1': ['I'],
        'chromosome_2': ['chr17'],
        'position_2': [7672484],
        'strand_2': ['-'],
        'operation_2': ['I'],
        'variant_size': [280],
        'variant_type': ['CRX'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-rna-008-tumor_ground_truth.tsv', sep='\t', index=False)

    # Step 6. Create transcript model alignments TSV file
    write_transcript_model_alignments_ground_truth(
        gencode=gencode,
        transcripts=[
            ('scga-mini-rna-008-tumor-1', tp53_tumor_structure, tp53_tumor_sequence),
            ('scga-mini-rna-008-tumor-2', tp53_normal_structure, tp53_normal_sequence)
        ],
        output_tsv_file='../../../../test/data/simulation/ground_truth/scga-mini-rna-008-tumor_transcript_model_alignments_ground_truth.tsv'
    )
