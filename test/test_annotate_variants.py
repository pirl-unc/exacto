from .data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import annotate_variant_calls


def test_annotate_variants_1():
    df_annotated = annotate_variant_calls(
        tsv_file=get_data_path(name='tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        gene_types=['protein_coding'],
        gene_levels=[1,2],
        transcript_types=['protein_coding'],
        transcript_levels=[1,2],
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert df_annotated['position_1_genic_region'].values[0] == 'exonic'
    assert df_annotated['position_2_genic_region'].values[0] == 'exonic'

