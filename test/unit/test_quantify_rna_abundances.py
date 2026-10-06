import pandas as pd
import pytest
from test.data import get_data_path
from exactolib.constants import OutputType
from exactolib.main import quantify_rna_abundances


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def _quantify_rna_abundances(sample_number: str, output_dir: str = '', output_type: OutputType = OutputType.DATAFRAME):
    # Pass-2 clusters; the consensus reference matches are call-rna-transcript-vars' on the stitched
    # consensus, whose transcripts keep their cluster IDs as names.
    prefix = 'assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_clusters' % sample_number
    return quantify_rna_abundances(
        clusters_tsv_file=get_data_path(name=prefix + '.tsv'),
        cluster_reference_matches_tsv_file=get_data_path(name=prefix + '_reference_transcripts.tsv'),
        cluster_splice_junctions_tsv_file=get_data_path(name=prefix + '_splice_junctions.tsv'),
        cluster_variants_tsv_file=get_data_path(name=prefix + '_variants_passed.tsv'),
        consensus_rna_reference_matches_tsv_file=get_data_path(name='variant_calling/rna/scga-mini-rna-%s-tumor_exacto_assembled_transcript_reference_transcript_matches.tsv' % sample_number),
        output_dir=output_dir,
        output_prefix='scga-mini-rna-%s-tumor' % sample_number,
        output_type=output_type
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_quantify_rna_abundances_counts_cluster_reads(sample_number):
    df_abundances, df_reference_transcript_abundances = _quantify_rna_abundances(sample_number=sample_number)
    prefix = 'assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_clusters' % sample_number
    df_clusters = pd.read_csv(get_data_path(name=prefix + '.tsv'), sep='\t')
    df_summary = pd.read_csv(get_data_path(name=prefix + '_summary.tsv'), sep='\t')
    df_variants = pd.read_csv(get_data_path(name=prefix + '_variants_passed.tsv'), sep='\t')

    # One row per cluster with its read count, shared reads included.
    df_abundances = df_abundances.sort_values('cluster_id')
    assert df_abundances['cluster_id'].tolist() == sorted(df_summary['cluster_id'])
    assert df_abundances['num_reads'].tolist() == df_summary.sort_values('cluster_id')['num_reads'].tolist()
    # The EM splits each read that sits in several clusters (e.g. one not covering the SNV of 001), so a
    # cluster's CPM lies between its exclusive reads and all of its reads over the N distinct reads,
    # and is n/N when it shares none.
    num_distinct_reads = df_clusters['read_name'].nunique()
    num_clusters_per_read = df_clusters.groupby('read_name')['cluster_id'].transform('nunique')
    for cluster_id, cpm in zip(df_abundances['cluster_id'], df_abundances['cpm']):
        is_cluster = df_clusters['cluster_id'] == cluster_id
        num_reads = int(is_cluster.sum())
        num_exclusive_reads = int((is_cluster & (num_clusters_per_read == 1)).sum())
        assert num_exclusive_reads / num_distinct_reads * 1e6 - 1e-6 <= cpm <= num_reads / num_distinct_reads * 1e6 + 1e-6
        if num_exclusive_reads == num_reads:
            assert cpm == pytest.approx(num_reads / num_distinct_reads * 1e6)
    assert df_abundances['cpm'].sum() == pytest.approx(1e6)
    # A cluster with a passed variant is the variant allele.
    assert set(df_abundances.loc[df_abundances['allele'] == 'variant', 'cluster_id']) == set(df_variants['cluster_id'])

    # The reference transcript table regroups the same reads.
    assert df_reference_transcript_abundances['num_clusters'].sum() == len(df_abundances)
    assert df_reference_transcript_abundances['num_reads'].sum() == df_abundances['num_reads'].sum()
    assert df_reference_transcript_abundances['cpm'].sum() == pytest.approx(1e6)


def test_quantify_rna_abundances_splits_tp53_by_allele():
    df_abundances, df_reference_transcript_abundances = _quantify_rna_abundances(sample_number='001')

    # scga-mini-rna-001: both clusters are TP53 ENST00000269305.9, cluster 1 with the SNV.
    df_abundances = df_abundances.sort_values('cluster_id')
    assert df_abundances[['cluster_id', 'reference_transcript_id', 'allele', 'num_reads']].values.tolist() == [
        [1, 'ENST00000269305.9', 'variant', 171],
        [2, 'ENST00000269305.9', 'reference', 174]
    ]
    assert df_reference_transcript_abundances[['reference_transcript_id', 'num_clusters', 'num_reads']].values.tolist() == [['ENST00000269305.9', 2, 345]]


def test_quantify_rna_abundances_writes_tables(tmp_path):
    _quantify_rna_abundances(sample_number='015', output_dir=str(tmp_path), output_type=OutputType.FILE)
    df_abundances, df_reference_transcript_abundances = _quantify_rna_abundances(sample_number='015')

    written = sorted(path.name for path in tmp_path.iterdir())
    assert len(written) == 2
    for name in written:
        df_file = pd.read_csv(tmp_path / name, sep='\t')
        df = df_abundances if 'cluster_id' in df_file.columns else df_reference_transcript_abundances
        assert list(df_file.columns) == list(df.columns), name
        assert len(df_file) == len(df), name
