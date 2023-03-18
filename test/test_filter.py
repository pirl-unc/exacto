import pandas as pd
from .data import get_data_path
from exacto.variant_filter import VariantFilter
from exacto.variants_list import VariantsList
from exacto.main import run_exacto_filter


def test_filter():
    # Step 1. Load data
    tsv_file = get_data_path(name='hg002_merged.tsv')
    germline_sv_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    hg38_gaps_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)
    df_excluded_variants = pd.read_csv(germline_sv_tsv_file, sep='\t')
    df_excluded_regions = pd.read_csv(hg38_gaps_tsv_file, sep='\t')

    # Step 2. Generate a list of VariantFilter instances
    variant_filter_1 = VariantFilter(
        quantifier='all',
        attribute='alt_tumor_reads',
        operator='>=',
        value=3
    )
    variant_filter_2 = VariantFilter(
        quantifier='all',
        attribute='filter',
        operator='==',
        value='"PASS"'
    )
    variant_filters = [variant_filter_1, variant_filter_2]

    # Step 2. Filter
    variants_list = run_exacto_filter(
        variants_list=variants_list,
        df_excluded_variants=df_excluded_variants,
        df_excluded_regions=df_excluded_regions,
        variant_filters=variant_filters,
        excluded_region_padding=100000,
        excluded_variant_padding=100,
        enforce_variant_type_checking=True
    )

    # Step 3. Write to file
    variants_list.to_dataframe().to_csv(
        get_data_path('hg002_merged_filtered.tsv'),
        sep='\t',
        index=False
    )
