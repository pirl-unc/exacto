import pandas as pd
import pytest

from test.data import get_data_path
from exactolib.main import identify_peptide_variants
from exactolib.variant_calling.peptide import _parse_intervals, identify_peptide_variants as identify_peptide_variants_


# scga-mini-rna-001: TP53 M246I, the one mutant amino acid at index 245 (counting from 0)
PROTEOFORMS_TSV_FILE = 'translation/scga-mini-rna-001-tumor_exacto_proteoforms.tsv'


def test_call_peptide_variants_reads_translate_transcripts_output():
    """The proteoforms TSV `translate-transcripts` writes is consumed as-is.

    Guards the column contract between the two commands: `proteoform_id` and
    `assembled_transcript_variant_ids` were renamed on the producing side, and
    reading the old names raised KeyError on every run.
    """
    df_peptide_variants = identify_peptide_variants(
        proteoforms_tsv_file=get_data_path(name=PROTEOFORMS_TSV_FILE),
        reference_proteome_fasta_file=get_data_path(name='exacto/exactolib/reference_peptides.fasta'),
        min_k=8,
        max_k=11,
        num_processes=1
    )
    df_proteoforms = pd.read_csv(get_data_path(name=PROTEOFORMS_TSV_FILE), sep='\t', keep_default_na=False)
    mutant_proteoform_id = df_proteoforms.loc[df_proteoforms['num_mutant_amino_acids'] > 0, 'proteoform_id'].tolist()

    # Every k from 8 to 11 has k windows over the one mutant position, none in the reference.
    assert len(df_peptide_variants) == 8 + 9 + 10 + 11
    # The output names the proteoform it came from, and carries its variant ids through.
    assert df_peptide_variants['proteoform_id'].unique().tolist() == mutant_proteoform_id
    assert 'assembled_transcript_variant_ids' in df_peptide_variants.columns
    for _, row in df_peptide_variants.iterrows():
        assert row['amino_acid_index_start'] <= 245 <= row['amino_acid_index_end']


def test_call_peptide_variants_drops_kmers_found_in_the_reference(tmp_path):
    # A reference holding the mutant TP53 stretch 236-255 removes every k-mer that lies inside it.
    df_proteoforms = pd.read_csv(get_data_path(name=PROTEOFORMS_TSV_FILE), sep='\t', keep_default_na=False)
    mutant_sequence = df_proteoforms.loc[df_proteoforms['num_mutant_amino_acids'] > 0, 'amino_acid_sequence'].iloc[0]
    reference_fasta_file = tmp_path / 'reference_peptides.fasta'
    reference_fasta_file.write_text('>ref_1\n%s\n' % mutant_sequence[236:256])

    df_peptide_variants = identify_peptide_variants(
        proteoforms_tsv_file=get_data_path(name=PROTEOFORMS_TSV_FILE),
        reference_proteome_fasta_file=str(reference_fasta_file),
        min_k=8,
        max_k=8,
        num_processes=1
    )

    # The eight 8-mers over 245 span 238-245 to 245-252, all inside 236-255.
    assert len(df_peptide_variants) == 0


def test_parse_intervals_accepts_both_encodings():
    """`format_intervals` emits `3:5;12;20:21`; older tables used `3-5,12,20-21`.

    Both must parse to the same positions. Reading the current encoding with the
    old parser raised ValueError on any proteoform with a range or more than one
    interval, so a single-interval fixture would not have caught it.
    """
    expected = {3, 4, 5, 12, 20, 21}
    assert _parse_intervals('3:5;12;20:21') == expected
    assert _parse_intervals('3-5,12,20-21') == expected

    assert _parse_intervals('245') == {245}
    assert _parse_intervals('') == set()
    assert _parse_intervals(None) == set()


def test_identify_peptide_variants_finds_kmers_across_a_multi_position_interval():
    """A mutant range yields every overlapping k-mer, not zero.

    Driven off a synthetic frame rather than the fixture: the real one carries a
    single mutant position, which cannot distinguish a working range parser from
    one that drops ranges on the floor.
    """
    df_proteoforms = pd.DataFrame({
        'proteoform_id': [7],
        'amino_acid_sequence': ['MKWVTFISLLFLFSSAYS'],
        'mutant_amino_acid_intervals': ['3:5'],
        'assembled_transcript_variant_ids': ['1;2'],
        'dna_variant_ids': ['']
    })

    df_peptide_variants = identify_peptide_variants_(
        df_proteoforms=df_proteoforms,
        reference_kmer_set={},
        min_k=3,
        max_k=3,
        num_processes=1
    )

    assert len(df_peptide_variants) > 0
    assert set(df_peptide_variants['proteoform_id']) == {7}
    assert set(df_peptide_variants['assembled_transcript_variant_ids']) == {'1;2'}
    # Every emitted 3-mer must overlap at least one of positions 3, 4, 5.
    for _, row in df_peptide_variants.iterrows():
        assert row['amino_acid_index_start'] <= 5 and row['amino_acid_index_end'] >= 3
