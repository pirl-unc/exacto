from .data import get_data_path
from exactolib.constants import GeneAnnotationSource, OutputType
from exactolib.main import integrate_variants


def test_integrate_variants_1():
    df_integrations = integrate_variants(
        dna_variant_call_annotation_set_tsv_file=get_data_path(name='tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv'),
        rna_variant_call_set_tsv_file=get_data_path(name='tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv'),
        reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
        reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
        reference_gene_annotation_assembly='hg38',
        reference_gene_annotation_version='v41',
        output_tsv_file='',
        output_type=OutputType.DATAFRAME
    )

    assert len(df_integrations) == 6
