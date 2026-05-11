import polars as pl
import pysam
import tempfile
from .data import get_data_path
from exactolib.main import build_genome_variation_graph
from exactolib.constants import GraphType


def test_build_genome_var_graph_1():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_dna_variant_callset_1.tsv')
    output_dir = tempfile.gettempdir()
    df_variants = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_dna_1.fasta'
    build_genome_variation_graph(
        df_variants=df_variants,
        fasta_file=get_data_path(name='fasta/sample.fa'),
        output_fasta_file=output_fasta_file,
        sequence_prefix='test_1',
        remove_unknown_bases=True,
        only_variant_sequences=False,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("test_1_1", 0, 20)
    sequence_2 = fasta.fetch("test_1_2", 0, 20)
    sequence_3 = fasta.fetch("test_1_3", 0, 20)

    assert sequence_1 == "ATGCATACGTAGCTAGCTAG"
    assert sequence_2 == "GGGTTTCCCAAAGGGTTTCC"
    assert sequence_3 == "GGATCGTATCTGACGTATGA"


def test_build_genome_var_graph_2():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_dna_variant_callset_2.tsv')
    output_dir = tempfile.gettempdir()
    df_variants = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_dna_2.fasta'
    build_genome_variation_graph(
        df_variants=df_variants,
        fasta_file=get_data_path(name='fasta/sample.fa'),
        output_fasta_file=output_fasta_file,
        sequence_prefix='test_2',
        remove_unknown_bases=True,
        only_variant_sequences=False,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("test_2_1", 0, 20)
    sequence_2 = fasta.fetch("test_2_2", 0, 20)
    sequence_3 = fasta.fetch("test_2_3", 0, 20)

    assert sequence_1 == "ATGCATACGTTAGCTAG"
    assert sequence_2 == "GGGTTTCCCAAAGGGTTTCC"
    assert sequence_3 == "GGATCGTATCTGACGTATGA"


def test_build_genome_var_graph_3():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_dna_variant_callset_3.tsv')
    output_dir = tempfile.gettempdir()
    df_variants = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_dna_3.fasta'
    build_genome_variation_graph(
        df_variants=df_variants,
        fasta_file=get_data_path(name='fasta/sample.fa'),
        output_fasta_file=output_fasta_file,
        sequence_prefix='test_3',
        remove_unknown_bases=True,
        only_variant_sequences=False,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("test_3_1", 0, 20)
    sequence_2 = fasta.fetch("test_3_2", 0, 20)
    sequence_3 = fasta.fetch("test_3_3", 0, 20)

    assert sequence_1 == "ATGCACGTACAGCTAGCTAG"
    assert sequence_2 == "GGGTTTCCCAAAGGGTTTCC"
    assert sequence_3 == "GGATCGTATCTGACGTATGA"


def test_build_genome_var_graph_4():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_dna_variant_callset_4.tsv')
    output_dir = tempfile.gettempdir()
    df_variants = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_dna_4.fasta'
    build_genome_variation_graph(
        df_variants=df_variants,
        fasta_file=get_data_path(name='fasta/sample.fa'),
        output_fasta_file=output_fasta_file,
        sequence_prefix='test_4',
        remove_unknown_bases=True,
        only_variant_sequences=False,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("test_4_1", 0, 10)
    sequence_2 = fasta.fetch("test_4_2", 0, 20)

    assert sequence_1 == "ATGCGTTTCC"
    assert sequence_2 == "GGATCGTATCTGACGTATGA"
