from .data import get_data_path
from exacto.variants_list import VariantsList
from exacto.main import run_exacto_merge


def test_merge():
    # Step 1. Load data
    tsv_file_1 = get_data_path(name='hg002_svim.tsv')
    tsv_file_2 = get_data_path(name='hg002_cutesv.tsv')
    tsv_file_3 = get_data_path(name='hg002_sniffles2.tsv')
    tsv_file_4 = get_data_path(name='hg002_pbsv.tsv')
    tsv_file_5 = get_data_path(name='hg002_deepvariant.tsv')

    variants_lists = []
    variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file_1))
    variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file_2))
    variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file_3))
    variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file_4))
    variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file_5))

    # Step 2. Merge
    variants_list = run_exacto_merge(
        variants_lists=variants_lists,
        enforce_variant_type_matching=True,
        max_neighbor_distance=1
    )

    # Step 3. Write to file
    df_variants = variants_list.to_dataframe()
    df_variants.to_csv(
        get_data_path('hg002_merged.tsv'),
        sep='\t',
        index=False
    )
