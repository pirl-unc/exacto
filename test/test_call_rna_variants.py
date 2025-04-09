import pandas as pd
from .data import get_data_path
from exactolib.main import identify_rna_variants


def test_call_rna_variants_1():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-100-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_2():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-101-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-101-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-101-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_3():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-102-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-102-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-102-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] == row['position_1']) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] == row['position_2']) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_4():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-103-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-103-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 2)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-103-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_5():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-104-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-104-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-104-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_6():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-105-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-105-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-105-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_7():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-106-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-106-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-106-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_8():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-107-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-107-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 2)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-107-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1


def test_call_rna_variants_9():
    df_exons,df_splice_junctions,df_variants = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-108-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-108-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        gene_annotation_source='gencode',
        output_exons_tsv_file='',
        output_sj_tsv_file='',
        output_variants_tsv_file='',
        gzip=False,
        output_type='dataframe'
    )

    assert(len(df_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-108-tumor_ground_truth.tsv'), sep='\t')
    for _,row in df_ground_truth.iterrows():
        df_matched = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        assert len(df_matched) == 1