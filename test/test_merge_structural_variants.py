from .data import get_data_path
from exactolib.main import *
from exactolib.utilities.merging_utils import *
from exactolib.constants import *
from exactolib.default_parameters import *


def test_merge_dna_structural_variants():
    # Step 1. Load data
    annotated_sniffles2_tsv_file = get_data_path(name='hg002_sniffles2_refined_annotated_ensembl.tsv')
    annotated_pbsv_tsv_file = get_data_path(name='hg002_pbsv_refined_annotated_ensembl.tsv')
    annotated_cutesv_tsv_file = get_data_path(name='hg002_cutesv_refined_annotated_ensembl.tsv')
    annotated_svim_tsv_file = get_data_path(name='hg002_svim_refined_annotated_ensembl.tsv')
    df_sniffles2 = pd.read_csv(annotated_sniffles2_tsv_file, sep='\t')
    df_pbsv = pd.read_csv(annotated_pbsv_tsv_file, sep='\t')
    df_cutesv = pd.read_csv(annotated_cutesv_tsv_file, sep='\t')
    df_svim = pd.read_csv(annotated_svim_tsv_file, sep='\t')

    # Step 2. Test
    list_df = [df_sniffles2, df_pbsv, df_svim, df_cutesv]
    df_merged, df_merged_deduped = run_exacto_merge_genomic_structural_variants(
        list_df=list_df,
        max_sv_cluster_distance=MAX_SV_CLUSTER_DISTANCE
    )

    # Step 3. Print output
    print("DataFrame columns:")
    print(df_merged_deduped.columns.values.tolist())
    print("DataFrame first row to dictionary:")
    print(df_merged_deduped.iloc[0].to_dict())

    # Step 4. Check for errors
    assert df_merged_deduped['variant_calling_methods'][0] == 'sniffles2,pbsv,svim,cutesv', \
        "The value of 'variant_calling_methods' of the first row is expected to be " \
        "'sniffles2,pbsv,svim,cutesv'. Instead it is '%s'." % (
        df_merged_deduped['variant_calling_methods'][0])

