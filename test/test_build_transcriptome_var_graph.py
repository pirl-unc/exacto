import polars as pl
import pysam
import tempfile
from .data import get_data_path
from exactolib.main import build_transcriptome_variation_graph
from exactolib.constants import GraphType


def test_build_transcriptome_var_graph_1():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_1.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_1.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, 20)

    assert sequence_1 == "ATGCGAGATAAGCGT"


def test_build_transcriptome_var_graph_2():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_2.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_2.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCCCCTAAGCGT"


def test_build_transcriptome_var_graph_3():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_3.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_3.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGATAAGCGT"


def test_build_transcriptome_var_graph_4():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_4.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_4.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTAGCAGCGT"


def test_build_transcriptome_var_graph_5():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_5.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_5.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTAAGCGTCGAGCACCAT"


def test_build_transcriptome_var_graph_6():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_6.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_6.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGCTAAGCGT"


def test_build_transcriptome_var_graph_7():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_7.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_7.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTCTCGCAT"


def test_build_transcriptome_var_graph_8():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_8.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_8.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTCTCGCCGGTCGA"


def test_build_transcriptome_var_graph_9():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_9.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_9.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAATCGC"


def test_build_transcriptome_var_graph_10():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_10.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_10.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "TAACTGCGATTACTG"


def test_build_transcriptome_var_graph_11():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_11.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_11.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "AGCTAAGCTA"


def test_build_transcriptome_var_graph_12():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_12.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_12.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTAAGCTAAGCGT"


def test_build_transcriptome_var_graph_13():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_13.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_13.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "TCGAAGCGTCGAGC"


def test_build_transcriptome_var_graph_14():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_14.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_14.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGGGGGAGCTAAGCGT"


def test_build_transcriptome_var_graph_15():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_15.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_15.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGCGCTAAGCGT"


def test_build_transcriptome_var_graph_16():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_16.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_16.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGTACGT"


def test_build_transcriptome_var_graph_17():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_17.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_17.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ACGTACGCAT"


def test_build_transcriptome_var_graph_18():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_18.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_18.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ATGCGAGCTAAGCGT"


def test_build_transcriptome_var_graph_19():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_19.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_19.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "ACGCTTAGCTCGCAT"


def test_build_transcriptome_var_graph_20():
    variant_tsv_file = get_data_path(name='tsv/variant_callset/sample_rna_variant_callset_20.tsv')
    fasta_file = get_data_path(name='fasta/sample3.fa')
    output_dir = tempfile.gettempdir()
    df_transcript_structures = pl.read_csv(variant_tsv_file, separator='\t')
    output_fasta_file = output_dir + '/test_rna_20.fasta'
    build_transcriptome_variation_graph(
        df_transcript_structures=df_transcript_structures,
        fasta_file=fasta_file,
        output_fasta_file=output_fasta_file,
        graph_type=str(GraphType.individual),
        num_threads=1
    )

    fasta = pysam.FastaFile(output_fasta_file)

    sequence_1 = fasta.fetch("1", 0, fasta.get_reference_length("1"))

    assert sequence_1 == "TAGCTCGCAT"

