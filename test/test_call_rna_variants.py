import pandas as pd
from .data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import identify_rna_variants


def filter_by_transcript_ids(df, expected_ids):
    expected = set(expected_ids.split(','))
    return df[df['reference_transcript_id'].apply(
        lambda s: set(s.split(',')) == expected
    )]

def test_call_rna_variants_1():
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_rna_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-100-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-101-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-101-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_rna_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-101-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-102-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-102-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_rna_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-102-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-103-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-103-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-103-tumor_ground_truth.tsv'), sep='\t')
    df_variants = filter_by_transcript_ids(
        df_rna_variants,
        'ENST00000698746.1,ENST00000570791.5'
    )
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-104-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-104-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_rna_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-104-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-105-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-105-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    assert(len(df_rna_variants) == 1)

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-105-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-106-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-106-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-106-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-107-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-107-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-107-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-108-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-108-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-108-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000269305.9']
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


def test_call_rna_variants_10():
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-109-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-109-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-109-tumor_ground_truth.tsv'), sep='\t')
    df_variants = filter_by_transcript_ids(
        df_rna_variants,
        'ENST00000333813.4,ENST00000570791.5,ENST00000263087.9'
    )
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


def test_call_rna_variants_11():
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-110-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-110-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-110-tumor_ground_truth.tsv'), sep='\t')
    df_variants = filter_by_transcript_ids(
        df_rna_variants,
        'ENST00000263092.11,ENST00000250113.12,ENST00000355530.7'
    )
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


def test_call_rna_variants_12():
    (df_assembled_transcripts,
     df_exons,
     df_introns,
     df_reference_transcript_matches,
     df_read_filter_status,
     df_transcript_model_structures,
     df_rna_variants) = identify_rna_variants(
        bam_file=get_data_path(name='bam/rna-111-tumor_minimap2_mdtagged_sorted.bam'),
        bam_bai_file=get_data_path(name='bam/rna-111-tumor_minimap2_mdtagged_sorted.bam.bai'),
        reference_genome_fasta_file=get_data_path(name='fasta/hg38_chr17-18.fa.gz'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        reference_transcript_top_k=3,
        output_dir='',
        output_prefix='',
        output_type=OutputType.DATAFRAME
    )

    df_ground_truth = pd.read_csv(get_data_path(name='tsv/ground_truth/rna-111-tumor_ground_truth.tsv'), sep='\t')
    df_variants = df_rna_variants[df_rna_variants['reference_transcript_id'] == 'ENST00000254719.10']
    for _,row in df_ground_truth.iterrows():
        df_matched_1 = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_1']) &
            (df_variants['position_1'] <= row['position_1'] + 100) &
            (df_variants['position_1'] >= row['position_1'] - 100) &
            (df_variants['operation_1'] == row['operation_1']) &
            (df_variants['chromosome_2'] == row['chromosome_2']) &
            (df_variants['position_2'] <= row['position_2'] + 100) &
            (df_variants['position_2'] >= row['position_2'] - 100) &
            (df_variants['operation_2'] == row['operation_2'])
        ]
        df_matched_2 = df_variants[
            (df_variants['chromosome_1'] == row['chromosome_2']) &
            (df_variants['position_1'] <= row['position_2'] + 100) &
            (df_variants['position_1'] >= row['position_2'] - 100) &
            (df_variants['operation_1'] == row['operation_2']) &
            (df_variants['chromosome_2'] == row['chromosome_1']) &
            (df_variants['position_2'] <= row['position_1'] + 100) &
            (df_variants['position_2'] >= row['position_1'] - 100) &
            (df_variants['operation_2'] == row['operation_1'])
        ]
        assert len(df_matched_1) == 1 or len(df_matched_2) == 1
