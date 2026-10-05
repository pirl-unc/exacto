import pandas as pd
import pytest
from test.data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import annotate_variant_calls


# (position_1_genic_region, position_2_genic_region) of each call, in table order, with the CLI's default
# filter (protein_coding genes and transcripts, levels 1 and 2). These are exacto-annotator's
# annotate_variant_calls_returns_matches_for_dna_NNN expectations except 006, whose Rust test runs
# unfiltered: there position 1 is exonic only in ENST00000698743.1, a nonsense_mediated_decay transcript.
EXPECTED_GENIC_REGIONS = {
    '001': [('exonic', 'exonic')],
    '002': [('exonic', 'exonic')],
    '003': [('exonic', 'exonic')],
    '004': [('exonic', 'exonic')],
    '005': [('exonic', 'exonic')],
    '006': [('intronic', 'intergenic'), ('intronic', 'intergenic')],
    '007': [('intronic', 'exonic')],
    '008': [('intronic', 'intronic')],
    '009': [('intronic', 'exonic')],
    '010': [('exonic', 'exonic')],
    '011': [('intronic', 'exonic')],
    '012': [('exonic', 'intronic')],
    '013': [('exonic', 'intronic'), ('intronic', 'intronic')],
    '014': [('exonic', 'exonic')],
    '015': [('intronic', 'exonic'), ('intronic', 'intronic')],
    '016': [('intronic', 'exonic')]
}

ANNOTATION_COLUMNS = [
    'position_1_genic_region', 'position_1_annotation',
    'position_1_plus_1_genic_region', 'position_1_plus_1_annotation',
    'position_2_minus_1_genic_region', 'position_2_minus_1_annotation',
    'position_2_genic_region', 'position_2_annotation'
]


@pytest.mark.parametrize('sample_number', sorted(EXPECTED_GENIC_REGIONS))
def test_annotate_variants_returns_genic_regions(sample_number):
    tsv_file = get_data_path(name='variant_calling/dna/scga-mini-dna-%s-tumor_exacto_germline_dna_variants.tsv' % sample_number)
    df_annotated = annotate_variant_calls(
        tsv_file=tsv_file,
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=['protein_coding'],
        gene_levels=[1,2],
        transcript_types=['protein_coding'],
        transcript_levels=[1,2],
        output_tsv_file='',
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )

    # The input columns and rows are kept, in order, with the annotation columns appended.
    df_variants = pd.read_csv(tsv_file, sep='\t')
    assert list(df_annotated.columns) == list(df_variants.columns) + ANNOTATION_COLUMNS
    assert df_annotated['variant_id'].tolist() == df_variants['variant_id'].tolist()

    genic_regions = list(zip(df_annotated['position_1_genic_region'], df_annotated['position_2_genic_region']))
    assert genic_regions == EXPECTED_GENIC_REGIONS[sample_number]


def test_annotate_variants_names_tp53_exon():
    # scga-mini-dna-001's SNV sits in TP53 (ENSG00000141510.18) exon ENSE00003712342.1 of ENST00000269305.9.
    df_annotated = annotate_variant_calls(
        tsv_file=get_data_path(name='variant_calling/dna/scga-mini-dna-001-tumor_exacto_germline_dna_variants.tsv'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=['protein_coding'],
        gene_levels=[1,2],
        transcript_types=['protein_coding'],
        transcript_levels=[1,2],
        output_tsv_file='',
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )

    assert 'ENSG00000141510.18' in df_annotated['position_1_annotation'].iloc[0]
    assert 'ENST00000269305.9' in df_annotated['position_1_annotation'].iloc[0]
    assert 'ENSE00003712342.1' in df_annotated['position_1_annotation'].iloc[0]


def test_annotate_variants_without_filters_includes_every_transcript_type():
    # Empty filter lists (the CLI flag given with no values) keep ENST00000698743.1, so 006 becomes exonic.
    df_annotated = annotate_variant_calls(
        tsv_file=get_data_path(name='variant_calling/dna/scga-mini-dna-006-tumor_exacto_germline_dna_variants.tsv'),
        reference_gene_annotation_file=get_data_path(name='references/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=[],
        gene_levels=[],
        transcript_types=[],
        transcript_levels=[],
        output_tsv_file='',
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )

    genic_regions = list(zip(df_annotated['position_1_genic_region'], df_annotated['position_2_genic_region']))
    assert genic_regions == [('exonic', 'intergenic'), ('exonic', 'intergenic')]
    assert 'ENST00000698743.1' in df_annotated['position_1_annotation'].iloc[0]
