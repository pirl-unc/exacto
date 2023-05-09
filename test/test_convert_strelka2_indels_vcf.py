from .data import get_data_path
from exacto.main import run_exacto_convert_vcf
from exacto.constants import *


def test_convert_strelka2_indels_vcf():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_hg001_strelka2_indels.vcf')

    # Step 2. Convert
    variants_list = run_exacto_convert_vcf(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.STRELKA2,
        sequencing_platform='pacbio'
    )

    # Step 3. Write to file
    df_variants = variants_list.to_dataframe()
    df_variants.to_csv(
        get_data_path('hg002_hg001_strelka2_indels.tsv'),
        sep='\t',
        index=False
    )