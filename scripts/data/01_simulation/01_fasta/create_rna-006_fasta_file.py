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

    # Step 3. Fetch WRAP53 sequence and create a somatic inversion
    df_transcript_wrap53 = gencode.df_transcripts[gencode.df_transcripts['transcript_id_stable'] == 'ENST00000698742']
    df_exons_wrap53 = gencode.df_exons[gencode.df_exons['transcript_id'] == df_transcript_wrap53['transcript_id'].values[0]]
    df_exons_wrap53.sort_values(by=['number'], ascending=False, inplace=True) # WRAP53 is on the forward strand

    # Step 4. Create an inverted WRAP53
    def get_exon_sequence(exon_number: int):
        df_matched = df_exons_wrap53[df_exons_wrap53['number'] == exon_number]
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

    intergenic_exon_1 = str(fasta.fetch('chr17', 7716000 - 1, 7716000 + 251))
    intergenic_exon_2 = str(fasta.fetch('chr17', 7715000 - 1, 7715000 + 251))
    intergenic_exon_3 = str(fasta.fetch('chr17', 7712000 - 1, 7712000 + 251))

    wrap53_tumor_sequence = \
        exon_1 + \
        exon_2 + \
        exon_3 + \
        exon_4 + \
        exon_5 + \
        reverse_complement(
            intergenic_exon_3 + \
            intergenic_exon_2 +
            intergenic_exon_1
        )

    wrap53_normal_sequence = \
        exon_1 + \
        exon_2 + \
        exon_3 + \
        exon_4 + \
        exon_5 + \
        exon_6 + \
        exon_7 + \
        exon_8 + \
        exon_9 + \
        exon_10 + \
        exon_11

    wrap53_tumor_peptide = longest_orf_peptide(seq=wrap53_tumor_sequence)
    wrap53_normal_peptide = longest_orf_peptide(seq=wrap53_normal_sequence)

    # The same transcripts as structures
    def get_exon_structure(exon_number: int):
        return TranscriptStructure.exon(fasta, df_exons_wrap53[df_exons_wrap53['number'] == exon_number].iloc[0])

    exon_structures = {number: get_exon_structure(exon_number=number) for number in range(1, 12)}
    intergenic_exon_1_structure = TranscriptStructure.segment(fasta, 'chr17', 7716000, 7716000 + 251)
    intergenic_exon_2_structure = TranscriptStructure.segment(fasta, 'chr17', 7715000, 7715000 + 251)
    intergenic_exon_3_structure = TranscriptStructure.segment(fasta, 'chr17', 7712000, 7712000 + 251)
    wrap53_tumor_structure = exon_structures[1] + exon_structures[2] + exon_structures[3] + exon_structures[4] + exon_structures[5] + \
        (intergenic_exon_3_structure + intergenic_exon_2_structure + intergenic_exon_1_structure).reverse_complement()
    wrap53_normal_structure = exon_structures[1] + exon_structures[2] + exon_structures[3] + exon_structures[4] + exon_structures[5] + exon_structures[6] + exon_structures[7] + exon_structures[8] + exon_structures[9] + exon_structures[10] + exon_structures[11]

    # Step 4. Create FASTA files
    create_fasta_file(
        sequences=[wrap53_tumor_sequence, wrap53_normal_sequence],
        sequence_names=['scga-mini-rna-006-tumor-1', 'scga-mini-rna-006-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-rna-006-tumor.fasta'
    )
    create_fasta_file(
        sequences=[wrap53_tumor_peptide, wrap53_normal_peptide],
        sequence_names=['scga-mini-pep-006-tumor-1', 'scga-mini-pep-006-tumor-2'],
        output_fasta_file='../../../../test/data/simulation/fasta/scga-mini-pep-006-tumor.fasta'
    )

    # Step 5. Create TSV file
    assert df_exons_wrap53.loc[df_exons_wrap53['number'] == 5, 'end'].values[0] == 7700829
    data = {
        'variant_id': [106],
        'chromosome_1': ['chr17'],
        'position_1': [7700829],
        'strand_1': ['+'],
        'operation_1': ['D'],
        'chromosome_2': ['chr17'],
        'position_2': [7716251],
        'strand_2': ['-'],
        'operation_2': ['D'],
        'variant_size': [''],
        'variant_type': ['BND'],
        'variant_sequence': ['']
    }
    pd.DataFrame(data).to_csv('../../../../test/data/simulation/ground_truth/scga-mini-rna-006-tumor_ground_truth.tsv', sep='\t', index=False)

    # Step 6. Create transcript model alignments TSV file
    write_transcript_model_alignments_ground_truth(
        gencode=gencode,
        transcripts=[
            ('scga-mini-rna-006-tumor-1', wrap53_tumor_structure, wrap53_tumor_sequence),
            ('scga-mini-rna-006-tumor-2', wrap53_normal_structure, wrap53_normal_sequence)
        ],
        output_tsv_file='../../../../test/data/simulation/ground_truth/scga-mini-rna-006-tumor_transcript_model_alignments_ground_truth.tsv'
    )
