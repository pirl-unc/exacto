import pandas as pd
import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import TranslationStrategy
from exactolib.main import translate_fastx_file, translate_sequence


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def test_translate_sequence_reads_rna_and_dna_alphabets():
    # The longest ORF starts at the first AUG (index 2) and ends after the UAA stop (index 13, exclusive).
    assert translate_sequence(rna_sequence='GGAUGAAAUGGUAAGG', strategy=TranslationStrategy.LONGEST_ORF) == [('MKW*', 2, 13)]
    assert translate_sequence(rna_sequence='GGATGAAATGGTAAGG', strategy=TranslationStrategy.LONGEST_ORF) == [('MKW*', 2, 13)]


def test_translate_sequence_without_start_codon_returns_nothing():
    assert translate_sequence(rna_sequence='GGCCCGGGUUU', strategy=TranslationStrategy.LONGEST_ORF) == []


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_translate_sequence_returns_simulated_peptides(sample_number):
    # scga-mini-pep-NNN-tumor-K is the longest ORF of scga-mini-rna-NNN-tumor-K, without its stop codon.
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-pep-%s-tumor.fasta' % sample_number)) as fasta:
        peptides = {record.name.replace('-pep-', '-rna-'): record.sequence for record in fasta}
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-rna-%s-tumor.fasta' % sample_number)) as fasta:
        for record in fasta:
            translations = translate_sequence(rna_sequence=record.sequence, strategy=TranslationStrategy.LONGEST_ORF)
            assert [peptide for peptide, _, _ in translations] == [peptides[record.name] + '*'], record.name


def test_translate_fastx_file_writes_peptides(tmp_path):
    # The FASTA is always written BGZF-compressed, whatever its name.
    output_fasta_file = str(tmp_path / 'scga-mini-rna-015-tumor_peptides.fasta.gz')
    output_tsv_file = str(tmp_path / 'scga-mini-rna-015-tumor_peptides.tsv')
    translate_fastx_file(
        fastx_file=get_data_path(name='simulation/fasta/scga-mini-rna-015-tumor.fasta'),
        output_fasta_file=output_fasta_file,
        output_tsv_file=output_tsv_file,
        strategy=TranslationStrategy.LONGEST_ORF,
        num_threads=1
    )

    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-pep-015-tumor.fasta')) as fasta:
        peptides = {record.name.replace('-pep-', '-rna-'): record.sequence + '*' for record in fasta}
    df_peptides = pd.read_csv(output_tsv_file, sep='\t')
    assert list(df_peptides.columns) == ['transcript_id', 'transcript_sequence', 'peptide_id', 'peptide_sequence', 'peptide_length', 'orf_start', 'orf_end']
    assert dict(zip(df_peptides['transcript_id'], df_peptides['peptide_sequence'])) == peptides
    with pysam.FastxFile(output_fasta_file) as fasta:
        assert [(record.name, record.sequence) for record in fasta] == list(zip(df_peptides['peptide_id'].astype(str), df_peptides['peptide_sequence']))
