from .data import get_data_path
from exacto.main import *
from exacto.utilities.merging_utils import *
from exacto.constants import *
from exacto.default_parameters import *


def test_merge_dna_small_variants():
    # Step 1. Load data
    refined_gatk4_mutect2_tsv_file = get_data_path(name='hg002_gatk4-mutect2_refined.tsv')
    refined_strelka2_tsv_file = get_data_path(name='hg002_strelka2_refined.tsv')

    df_mutect2 = pd.read_csv(refined_gatk4_mutect2_tsv_file, sep='\t')
    df_strelka2 = pd.read_csv(refined_strelka2_tsv_file, sep='\t')

    # Step 2. Test
    list_df = [df_mutect2, df_strelka2]
    df_merged, df_merged_deduped = run_exacto_merge_genomic_small_variants(
        list_df=list_df,
        enforce_variant_type_matching=True,
        max_clustering_distance=MAX_SMALL_VARIANT_CLUSTER_DISTANCE
    )

    # Step 3. Print output
    print("DataFrame columns:")
    print(df_merged_deduped.columns.values.tolist())
    print("DataFrame first row to dictionary:")
    print(df_merged_deduped.iloc[0].to_dict())

    # Step 4. Write to file
    output_tsv_file = get_data_path('hg002_small_variants_merged_deduped.tsv')
    df_merged_deduped.to_csv(output_tsv_file, sep='\t', index=False)
