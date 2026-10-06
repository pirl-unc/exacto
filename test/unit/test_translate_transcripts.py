import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import TranslationStrategy, OutputType
from exactolib.main import translate_transcripts


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def _translate_transcripts(sample_number: str, **overrides):
    # The scga-mini inputs scripts/data/07_translation/01_run_exacto_translate_transcripts.sh passes
    arguments = dict(
        assembled_transcript_model_alignments_tsv_file=get_data_path(name='variant_calling/rna/scga-mini-rna-%s-tumor_exacto_assembled_transcript_model_alignments.tsv' % sample_number),
        assembled_transcript_variants_tsv_file=get_data_path(name='variant_calling/rna/scga-mini-rna-%s-tumor_exacto_assembled_transcript_variants.tsv' % sample_number),
        dna_variants_tsv_file=get_data_path(name='variant_calling/dna/scga-mini-dna-%s-tumor_exacto_somatic_dna_variants.tsv' % sample_number),
        integrated_variants_tsv_file=get_data_path(name='integration/scga-mini-rna-%s-tumor_exacto_integrated_variants.tsv' % sample_number),
        rna_consensus_tsv_file=get_data_path(name='assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus.tsv' % sample_number),
        stitched_transcripts_tsv_file=get_data_path(name='assembly/rna/stitched/scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts.tsv' % sample_number),
        strategy=TranslationStrategy.LONGEST_ORF,
        output_dir='',
        output_prefix='',
        num_threads=2,
        output_type=OutputType.DATAFRAME
    )
    arguments.update(overrides)
    return translate_transcripts(**arguments)


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_translate_transcripts_returns_simulated_peptides(sample_number):
    df_proteoforms, df_proteoform_nucleotides = _translate_transcripts(sample_number=sample_number)

    # One proteoform per simulated transcript, each the simulated peptide (stop codons aside).
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-pep-%s-tumor.fasta' % sample_number)) as fasta:
        expected_sequences = sorted(record.sequence.replace('*', '') for record in fasta)
    assert sorted(df_proteoforms['amino_acid_sequence'].str.replace('*', '', regex=False)) == expected_sequences
    assert set(df_proteoform_nucleotides['proteoform_id']) == set(df_proteoforms['proteoform_id'])


def test_translate_transcripts_finds_tp53_m246i():
    df_proteoforms, df_proteoform_nucleotides = _translate_transcripts(sample_number='001')

    # scga-mini-rna-001: TP53 M246I (MGGMN -> MGGIN), amino acid 245 counting from 0, from DNA variant 1.
    df_mutant = df_proteoforms[df_proteoforms['num_mutant_amino_acids'] > 0]
    assert len(df_mutant) == 1
    assert df_mutant['reference_gene_name'].iloc[0] == 'TP53'
    assert df_mutant['mutant_amino_acid_intervals'].iloc[0] == '245:245'
    assert df_mutant['amino_acid_sequence'].iloc[0][242:247] == 'MGGIN'
    assert str(df_mutant['dna_variant_ids'].iloc[0]) == '1'

    df_nucleotides = df_proteoform_nucleotides[
        (df_proteoform_nucleotides['proteoform_id'] == df_mutant['proteoform_id'].iloc[0]) &
        df_proteoform_nucleotides['is_amino_acid_variant']
    ]
    assert set(df_nucleotides['amino_acid_index']) == {245}
    assert set(df_nucleotides['amino_acid']) == {'I'}


def test_translate_transcripts_writes_tables(tmp_path):
    df_proteoforms, df_proteoform_nucleotides = _translate_transcripts(sample_number='001', output_dir=str(tmp_path), output_prefix='scga-mini-rna-001-tumor', output_type=OutputType.FILE)

    # File mode writes the tables and the proteoform FASTA, and returns empty frames.
    assert len(df_proteoforms) == 0 and len(df_proteoform_nucleotides) == 0
    for name in ['scga-mini-rna-001-tumor_exacto_proteoforms.tsv',
                 'scga-mini-rna-001-tumor_exacto_proteoform_nucleotides.tsv',
                 'scga-mini-rna-001-tumor_exacto_proteoforms.fasta']:
        assert (tmp_path / name).exists(), name


def test_translate_transcripts_requires_exactly_one_sequence_source():
    """Transcript sequences and their read names must come from a single file.

    Supplying neither leaves the translator with no sequences at all; supplying both invites
    two files that disagree about what a transcript is. Checked in `main` rather than only in
    argparse so the Python API is guarded for callers that never touch the CLI.
    """
    consensus_tsv_file = get_data_path(name='assembly/rna/pass2/scga-mini-rna-001-tumor_pass2_exacto_rna_consensus.tsv')

    with pytest.raises(ValueError):
        _translate_transcripts(sample_number='001', rna_consensus_tsv_file='')

    with pytest.raises(ValueError):
        _translate_transcripts(sample_number='001', assembled_transcript_support_tsv_file=consensus_tsv_file)
