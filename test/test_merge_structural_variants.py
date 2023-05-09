from .data import get_data_path
from exacto.main import run_exacto_convert_vcf, run_exacto_merge_variant_calls
from exacto.constants import VariantCallingMethods
from exacto.default_parameters import MERGE_MAX_NEIGHBOR_DISTANCE


def test_merge_structural_variants():
    # Step 1. Load data
    vcf_file_1 = get_data_path(name='hg002_cutesv.vcf')
    vcf_file_2 = get_data_path(name='hg002_pbsv.vcf')
    vcf_file_3 = get_data_path(name='hg002_sniffles2.vcf')
    vcf_file_4 = get_data_path(name='hg002_svim.vcf')

    # Step 2. Convert
    variants_list_1 = run_exacto_convert_vcf(
        vcf_file=vcf_file_1,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.CUTESV,
        sequencing_platform='pacbio'
    )
    variants_list_2 = run_exacto_convert_vcf(
        vcf_file=vcf_file_2,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.PBSV,
        sequencing_platform='pacbio'
    )
    variants_list_3 = run_exacto_convert_vcf(
        vcf_file=vcf_file_3,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SNIFFLES2,
        sequencing_platform='pacbio'
    )
    variants_list_4 = run_exacto_convert_vcf(
        vcf_file=vcf_file_4,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SVIM,
        sequencing_platform='pacbio'
    )

    # Step 3. Merge
    variants_lists = [
        variants_list_1,
        variants_list_2,
        variants_list_3,
        variants_list_4
    ]
    variants_list_merged = run_exacto_merge_variant_calls(
        variants_lists=variants_lists,
        max_neighbor_distance=MERGE_MAX_NEIGHBOR_DISTANCE
    )

    # Step 4. Write to file
    df_variants = variants_list_merged.to_dataframe()
    df_variants.to_csv(
        get_data_path('hg002_structural_variants_merged.tsv'),
        sep='\t',
        index=False
    )
