from .data import get_data_path
from exactolib.main import identify_peptide_variants


def test_call_peptide_variants_1():
    df_pep_variants = identify_peptide_variants(
        primary_structures_tsv_file=get_data_path(name="tsv/primary_structure/rna-100-tumor_minimap2_mdtagged_sorted_exacto_primary_structures.tsv"),
        reference_proteome_fasta_file=get_data_path(name='fasta/reference_peptides.fasta'),
        min_k=3,
        max_k=8,
        num_processes=1
    )
    assert len(df_pep_variants) != 0
