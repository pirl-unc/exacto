from .data import get_data_path
from exactolib.constants import TranslationStrategy, OutputType
from exactolib.main import translate_structures


def test_translate_structs_1():
    df_primary_structures = translate_structures(
        transcript_structures_tsv_file=get_data_path(name='tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv'),
        rna_variant_calls_tsv_file=get_data_path(name='tsv/variant_callset/rna-100-tumor-minimap2_mdtagged_sorted_bam_exacto_call_rna_vars_outputs/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv'),
        integrated_variants_tsv_file=get_data_path(name='tsv/variant_callset/rna-100-tumor_dna-001-tumor_variants_integrated.tsv'),
        strategy=TranslationStrategy.LONGEST_ORF,
        output_tsv_file='',
        output_fasta_file='',
        num_threads=1,
        output_type=OutputType.DATAFRAME,
    )

    print(df_primary_structures.head(n=10000))

    # assert(len(df_peptides) == 2)
