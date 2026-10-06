import pandas as pd
import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import OutputType
from exactolib.main import determine_rna_consensus


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def _determine_rna_consensus(sample_number: str, output_dir: str = '', output_type: OutputType = OutputType.DATAFRAME):
    # As scripts/data/04_assembly/rna/05_run_exacto_determine_rna_consensus.sh runs it: pass-2 clusters
    # of the pass-1 corrected reads
    return determine_rna_consensus(
        tsv_file=get_data_path(name='assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_clusters.tsv' % sample_number),
        fastq_file=get_data_path(name='assembly/rna/pass1/scga-mini-rna-%s-tumor_pass1_exacto_rna_corrected_reads.fastq.gz' % sample_number),
        output_dir=output_dir,
        output_prefix='scga-mini-rna-%s-tumor_pass2' % sample_number,
        num_threads=2,
        output_type=output_type
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_determine_rna_consensus_returns_simulated_transcripts(sample_number):
    df_consensus = _determine_rna_consensus(sample_number=sample_number)

    # One consensus per cluster, built from all of its reads.
    df_summary = pd.read_csv(get_data_path(name='assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_clusters_summary.tsv' % sample_number), sep='\t')
    assert sorted(df_consensus.select(['cluster_id', 'num_reads']).rows()) == sorted(df_summary[['cluster_id', 'num_reads']].itertuples(index=False, name=None))
    for read_names, num_reads in df_consensus.select(['read_names', 'num_reads']).rows():
        assert len(read_names.split(';')) == num_reads

    # Each consensus is the full-length simulated transcript of its cluster.
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-rna-%s-tumor.fasta' % sample_number)) as fasta:
        simulated_sequences = {record.sequence.upper(): record.name for record in fasta}
    names = []
    for sequence in df_consensus['consensus_sequence'].str.to_uppercase().to_list():
        name = simulated_sequences.get(sequence, simulated_sequences.get(_reverse_complement(sequence)))
        assert name is not None
        names.append(name)
    assert sorted(names) == sorted(simulated_sequences.values())


def test_determine_rna_consensus_writes_tsv_and_fasta(tmp_path):
    _determine_rna_consensus(sample_number='015', output_dir=str(tmp_path), output_type=OutputType.FILE)
    df_consensus = _determine_rna_consensus(sample_number='015')

    # The FASTA (what minimap2 realigns before stitch-reference-transcripts) names each consensus by its cluster.
    df_file = pd.read_csv(tmp_path / 'scga-mini-rna-015-tumor_pass2_exacto_rna_consensus.tsv', sep='\t')
    assert list(df_file.columns) == df_consensus.columns
    assert df_file['consensus_sequence'].tolist() == df_consensus['consensus_sequence'].to_list()
    with pysam.FastxFile(str(tmp_path / 'scga-mini-rna-015-tumor_pass2_exacto_rna_consensus.fasta')) as fasta:
        assert [(record.name, record.sequence) for record in fasta] == list(zip(df_file['cluster_id'].astype(str), df_file['consensus_sequence']))
