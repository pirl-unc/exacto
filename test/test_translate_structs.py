from .data import get_data_path
from exactolib.constants import TranslationStrategy, OutputType
from exactolib.main import translate_structures


def test_translate_structs_1():
    df_primary_structures = translate_structures(
        rna_assembly_support_tsv_file=get_data_path(name='tsv/rna_assembly/rna-100_transcriptome_assembly_read_support.tsv'),
        transcript_model_structures_tsv_file=get_data_path(name='tsv/transcript_model_structure/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv'),
        dna_variants_tsv_file=get_data_path(name='tsv/variant_callset/integrate_variants/dna-001-tumor_minimap2_mdtagged_sorted_exacto_dna_variant_records.tsv'),
        rna_variants_tsv_file=get_data_path(name='tsv/variant_callset/integrate_variants/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_records.tsv'),
        integrated_variants_tsv_file=get_data_path(name='tsv/variant_callset/rna-100_dna-001_integration.tsv'),
        strategy=TranslationStrategy.LONGEST_ORF,
        output_dir='',
        output_prefix='',
        num_threads=1,
        output_type=OutputType.DATAFRAME
    )
