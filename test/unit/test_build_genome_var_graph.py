import polars as pl
import pysam
import tempfile
from pathlib import Path
from test.data import get_data_path
from exactolib.main import build_genome_variation_graph
from exactolib.constants import GraphType, OutputType


# sample.fa holds chrA, chrB and chrC; with only_variant_sequences=False every contig's reference
# sequence follows the variant sequences, in FASTA order.
REFERENCE_SEQUENCES = [
    "ATGCGTACGTAGCTAGCTAG",
    "GGGTTTCCCAAAGGGTTTCC",
    "GGATCGTATCTGACGTATGA"
]


def _build_records(callset: int, only_variant_sequences: bool = False):
    with tempfile.TemporaryDirectory() as temp_dir:
        output_fasta_file = str(Path(temp_dir) / ('test_dna_%i.fasta' % callset))
        build_genome_variation_graph(
            df_variants=pl.read_csv(get_data_path(name='exacto/exacto-graph/sample_dna_variant_callset_%i.tsv' % callset), separator='\t'),
            fasta_file=get_data_path(name='exacto/exacto-graph/sample.fa'),
            output_fasta_file=output_fasta_file,
            sequence_prefix='test_%i' % callset,
            remove_unknown_bases=True,
            only_variant_sequences=only_variant_sequences,
            graph_type=str(GraphType.individual),
            num_threads=1
        )
        with pysam.FastxFile(output_fasta_file) as fasta:
            return [(record.name, record.sequence) for record in fasta]


def test_build_genome_var_graph_1():
    records = _build_records(callset=1)

    assert [name for name, _ in records] == ["test_1_1", "test_1_2", "test_1_3", "test_1_4"]
    assert [sequence for _, sequence in records] == ["ATGCATACGTAGCTAGCTAG"] + REFERENCE_SEQUENCES


def test_build_genome_var_graph_2():
    records = _build_records(callset=2)

    assert [sequence for _, sequence in records] == ["ATGCATACGTTAGCTAG"] + REFERENCE_SEQUENCES


def test_build_genome_var_graph_3():
    records = _build_records(callset=3)

    assert [sequence for _, sequence in records] == ["ATGCACGTACAGCTAGCTAG"] + REFERENCE_SEQUENCES


def test_build_genome_var_graph_4():
    records = _build_records(callset=4)

    assert [sequence for _, sequence in records] == ["ATGCGTTTCC"] + REFERENCE_SEQUENCES


def test_build_genome_var_graph_only_variant_sequences():
    records = _build_records(callset=1, only_variant_sequences=True)

    assert records == [("test_1_1", "ATGCATACGTAGCTAGCTAG")]


def test_build_genome_var_graph_dataframe_flags_variant_sequences():
    df_sequences = build_genome_variation_graph(
        df_variants=pl.read_csv(get_data_path(name='exacto/exacto-graph/sample_dna_variant_callset_4.tsv'), separator='\t'),
        fasta_file=get_data_path(name='exacto/exacto-graph/sample.fa'),
        output_fasta_file='',
        sequence_prefix='test_4',
        remove_unknown_bases=True,
        only_variant_sequences=False,
        graph_type=str(GraphType.individual),
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )

    assert list(df_sequences.columns) == ['id', 'sequence', 'is_variant']
    assert df_sequences['sequence'].tolist() == ["ATGCGTTTCC"] + REFERENCE_SEQUENCES
    assert df_sequences['is_variant'].tolist() == [True, False, False, False]
