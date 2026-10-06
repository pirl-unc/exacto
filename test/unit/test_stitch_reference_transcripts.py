import pandas as pd
import pysam
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import stitch_reference_transcripts


SAMPLE_NUMBERS = ['%03i' % i for i in range(1, 17)]


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def _stitch_reference_transcripts(sample_number: str, output_dir: str = '', output_type: OutputType = OutputType.DATAFRAME) -> pd.DataFrame:
    # The realigned pass-2 consensus, as scripts/data/04_assembly/rna/07_run_exacto_stitch_reference_transcripts.sh runs it
    return stitch_reference_transcripts(
        bam_file=get_data_path(name='assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus_minimap2_sorted.bam' % sample_number),
        reference_genome_fasta_file=get_data_path(name='references/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir=output_dir,
        output_prefix='scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus' % sample_number,
        num_threads=2,
        output_type=output_type
    )


@pytest.mark.parametrize('sample_number', SAMPLE_NUMBERS)
def test_stitch_reference_transcripts_returns_simulated_transcripts(sample_number):
    df_stitched = _stitch_reference_transcripts(sample_number=sample_number)

    # One row per consensus sequence, each the full-length simulated transcript it came from.
    with pysam.FastxFile(get_data_path(name='simulation/fasta/scga-mini-rna-%s-tumor.fasta' % sample_number)) as fasta:
        simulated_sequences = {record.sequence.upper(): record.name for record in fasta}
    with pysam.FastxFile(get_data_path(name='assembly/rna/pass2/scga-mini-rna-%s-tumor_pass2_exacto_rna_consensus.fasta' % sample_number)) as fasta:
        consensus_names = sorted(record.name for record in fasta)

    assert sorted(df_stitched['read_name'].astype(str)) == consensus_names
    names = []
    for sequence in df_stitched['stitched_sequence'].str.upper():
        name = simulated_sequences.get(sequence, simulated_sequences.get(_reverse_complement(sequence)))
        assert name is not None
        names.append(name)
    assert sorted(names) == sorted(simulated_sequences.values())
    assert (df_stitched['stitched_length'] == df_stitched['stitched_sequence'].str.len()).all()


def test_stitch_reference_transcripts_writes_tsv_and_fasta(tmp_path):
    df_empty = _stitch_reference_transcripts(sample_number='007', output_dir=str(tmp_path), output_type=OutputType.FILE)
    df_stitched = _stitch_reference_transcripts(sample_number='007')

    # File mode writes the table and the FASTA translate-transcripts and minimap2 read, and returns nothing.
    assert len(df_empty) == 0
    prefix = 'scga-mini-rna-007-tumor_pass2_exacto_rna_consensus_exacto_stitched_transcripts'
    df_file = pd.read_csv(tmp_path / ('%s.tsv' % prefix), sep='\t', dtype={'read_name': str})
    assert list(df_file.columns) == list(df_stitched.columns)
    assert df_file['stitched_sequence'].tolist() == df_stitched['stitched_sequence'].tolist()
    with pysam.FastxFile(str(tmp_path / ('%s.fasta' % prefix))) as fasta:
        assert [(record.name, record.sequence) for record in fasta] == list(zip(df_file['read_name'], df_file['stitched_sequence']))
