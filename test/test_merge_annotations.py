from .data import get_data_path
from exacto.main import *
from exacto.utilities.merging_utils import *
from exacto.constants import *
from exacto.default_parameters import *


def test_merge_annotations():
    # Step 1. Load data
    gatk4_mutect2_gencode_tsv_file = get_data_path(name='hg002_gatk4-mutect2_annotated_gencode.tsv')
    gatk4_mutect2_ensembl_tsv_file = get_data_path(name='hg002_gatk4-mutect2_annotated_ensembl.tsv')

    df_mutect2_gencode = pd.read_csv(gatk4_mutect2_gencode_tsv_file, sep='\t')
    df_mutect2_ensembl = pd.read_csv(gatk4_mutect2_ensembl_tsv_file, sep='\t')

    # Step 2. Test
    list_df = [df_mutect2_gencode, df_mutect2_ensembl]
    df_merged = run_exacto_merge_annotations(list_df=list_df)

    # Step 4. Write to file
    output_tsv_file = get_data_path('hg002_small_variants_gencode_ensembl_merged_deduped.tsv')
    df_merged.to_csv(output_tsv_file, sep='\t', index=False)
