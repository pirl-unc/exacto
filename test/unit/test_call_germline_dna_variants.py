import pandas as pd
import pytest
from test.data import get_data_path
from exactolib.constants import OutputType
from exactolib.default import CALL_DNA_VARS_PRESETS
from exactolib.main import identify_germline_dna_variants


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]

# The two copies of the 012 duplication meet across a G that both sides hold, so the junction can be
# written one base either way (the caller writes the G as inserted sequence). Same rule as exacto-caller.
POSITION_TOLERANCE = {'012': 1}


def _identify_germline_dna_variants(sample_number: str, regions=[]) -> pd.DataFrame:
    # The tumor BAM alone with --preset pb, as scripts/data/05_variant_calling/01_run_exacto_call_germline_dna_vars.sh runs it
    return identify_germline_dna_variants(
        bam_file=get_data_path(name='alignment/scga-mini-dna-%s-tumor_minimap2_sorted.bam' % sample_number),
        bam_bai_file=get_data_path(name='alignment/scga-mini-dna-%s-tumor_minimap2_sorted.bam.bai' % sample_number),
        fasta_file=get_data_path(name='references/hg38_chr17-18.fa.gz'),
        output_tsv_file='',
        regions=regions,
        num_threads=4,
        output_type=OutputType.DATAFRAME,
        **CALL_DNA_VARS_PRESETS['pb']
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_identify_germline_dna_variants_matches_ground_truth(sample_number):
    df_variants = _identify_germline_dna_variants(sample_number=sample_number)
    df_ground_truth = pd.read_csv(get_data_path(name='simulation/ground_truth/scga-mini-dna-%s-tumor_ground_truth.tsv' % sample_number), sep='\t', keep_default_na=False)
    tolerance = POSITION_TOLERANCE.get(sample_number, 0)

    # Each ground truth row is exactly one call, and there are no other calls.
    for _, row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            df_variants['position_1'].between(row['position_1'] - tolerance, row['position_1'] + tolerance) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            df_variants['position_2'].between(row['position_2'] - tolerance, row['position_2'] + tolerance) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1, 'ground truth row %i' % row['variant_call_id']
        if tolerance == 0:
            assert df_matched['sequence'].iloc[0].upper() == row['variant_sequence'].upper()
    assert len(df_variants) == len(df_ground_truth)


def test_identify_germline_dna_variants_dataframe_matches_file_output(tmp_path):
    output_tsv_file = str(tmp_path / 'scga-mini-dna-001-tumor_exacto_germline_dna_variants.tsv')
    identify_germline_dna_variants(
        bam_file=get_data_path(name='alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam'),
        bam_bai_file=get_data_path(name='alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai'),
        fasta_file=get_data_path(name='references/hg38_chr17-18.fa.gz'),
        output_tsv_file=output_tsv_file,
        regions=[],
        num_threads=4,
        output_type=OutputType.FILE,
        **CALL_DNA_VARS_PRESETS['pb']
    )
    df_file = pd.read_csv(output_tsv_file, sep='\t', keep_default_na=False)
    df_variants = _identify_germline_dna_variants(sample_number='001').astype(str)

    assert list(df_file.columns) == list(df_variants.columns)
    assert df_file.astype(str).values.tolist() == df_variants.values.tolist()


def test_identify_germline_dna_variants_restricted_to_a_region():
    # scga-mini-dna-001's SNV is chr17:7674225 (TP53); a region elsewhere on chr17 holds no call.
    assert len(_identify_germline_dna_variants(sample_number='001', regions=[('chr17', 7_670_000, 7_680_000)])) == 1
    assert len(_identify_germline_dna_variants(sample_number='001', regions=[('chr17', 1_000_000, 2_000_000)])) == 0
