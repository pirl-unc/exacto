from .data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import integrate_variants


def test_integrate_variants_1():
    df_integrations = integrate_variants(
        dna_variants_tsv_file=get_data_path(name='tsv/variant_callset/integrate_variants/dna-001-tumor_minimap2_mdtagged_sorted_exacto_dna_variant_records.tsv'),
        rna_variants_tsv_file=get_data_path(name='tsv/variant_callset/integrate_variants/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_records.tsv'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert len(df_integrations) == 6
