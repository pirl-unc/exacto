import pandas as pd
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import integrate_variants


# scga-mini-dna-014 has no somatic call (see test_call_somatic_dna_variants.py), so nothing integrates.
SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17) if i != 14]

# The RNA variant is the DNA variant itself (or, for 005, the splice site it deletes).
SAME_SITE_SAMPLE_NUMBERS = ['001', '002', '003', '004', '005', '010']

COLUMNS = [
    'assembled_transcript_name', 'reference_gene_name', 'reference_transcript_id', 'rna_variant_id',
    'dna_variant_id', 'distance', 'rna_variant_position', 'dna_variant_position'
]


def _integrate_variants(sample_number: str) -> pd.DataFrame:
    # As scripts/data/06_integration/01_run_exacto_integrate_vars.sh runs it: the DNA sample of the same number
    return integrate_variants(
        dna_variants_tsv_file=get_data_path(name='variant_calling/dna/scga-mini-dna-%s-tumor_exacto_somatic_dna_variants.tsv' % sample_number),
        rna_variants_tsv_file=get_data_path(name='variant_calling/rna/scga-mini-rna-%s-tumor_exacto_assembled_transcript_variants.tsv' % sample_number),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_tsv_file='',
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_integrate_variants_links_every_dna_variant(sample_number):
    df_integrations = _integrate_variants(sample_number=sample_number)
    df_dna_variants = pd.read_csv(get_data_path(name='variant_calling/dna/scga-mini-dna-%s-tumor_exacto_somatic_dna_variants.tsv' % sample_number), sep='\t')
    df_rna_variants = pd.read_csv(get_data_path(name='variant_calling/rna/scga-mini-rna-%s-tumor_exacto_assembled_transcript_variants.tsv' % sample_number), sep='\t')

    assert list(df_integrations.columns) == COLUMNS
    # Each simulated DNA event reaches the transcripts made from it, and every row names real calls.
    assert set(df_integrations['dna_variant_id']) == set(df_dna_variants['variant_id'])
    assert set(df_integrations['rna_variant_id']) <= set(df_rna_variants['variant_id'])


@pytest.mark.parametrize('sample_number', SAME_SITE_SAMPLE_NUMBERS)
def test_integrate_variants_links_same_site_variants_at_distance_zero(sample_number):
    df_integrations = _integrate_variants(sample_number=sample_number)

    assert df_integrations[['rna_variant_id', 'dna_variant_id', 'distance']].values.tolist() == [[1, 1, 0]]
    assert df_integrations[['rna_variant_position', 'dna_variant_position']].values.tolist() == [['position_1', 'position_1']]


def test_integrate_variants_returns_header_for_dna_table_without_rows():
    df_integrations = _integrate_variants(sample_number='014')

    assert list(df_integrations.columns) == COLUMNS
    assert len(df_integrations) == 0


def test_integrate_variants_writes_dataframe_to_file(tmp_path):
    output_tsv_file = str(tmp_path / 'scga-mini-rna-006-tumor_exacto_integrated_variants.tsv')
    integrate_variants(
        dna_variants_tsv_file=get_data_path(name='variant_calling/dna/scga-mini-dna-006-tumor_exacto_somatic_dna_variants.tsv'),
        rna_variants_tsv_file=get_data_path(name='variant_calling/rna/scga-mini-rna-006-tumor_exacto_assembled_transcript_variants.tsv'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_tsv_file=output_tsv_file,
        num_threads=1,
        output_type=OutputType.FILE
    )
    df_file = pd.read_csv(output_tsv_file, sep='\t', dtype={'assembled_transcript_name': str})
    df_integrations = _integrate_variants(sample_number='006')

    assert list(df_file.columns) == COLUMNS
    assert df_file.values.tolist() == df_integrations.values.tolist()
