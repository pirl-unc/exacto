# from .data import get_data_path
# from exactolib.constants import GeneAnnotationSource, OutputType
# from exactolib.main import integrate_dna_rna_variants
#
#
# def test_integrate_variants_1():
#     df_integrations = integrate_dna_rna_variants(
#         annotated_dna_variant_callset_tsv_file=get_data_path(name='tsv/variant_callset/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv'),
#         rna_variant_callset_tsv_file=get_data_path(name='tsv/variant_callset/rna-100-tumor_minimap2_mdtagged_sorted_exacto_variants.tsv'),
#         reference_gene_annotation_file=get_data_path(name='gtf/gencode.v41.annotation.chr17-18.gtf.gz'),
#         reference_gene_annotation_source=GeneAnnotationSource.GENCODE,
#         output_tsv_file='',
#         output_type=OutputType.DATAFRAME
#     )
#
#     assert len(df_integrations) == 4
