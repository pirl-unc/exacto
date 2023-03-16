from .data import get_data_path
from exacto.main import run_exacto_convert
from exacto.constants import *


def test_convert_deepvariant_vcf():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_deepvariant.vcf')

    # Step 2. Convert
    variants_list = run_exacto_convert(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.DEEPVARIANT,
        sequencing_platform='pacbio',
        tumor_sample_id='hg002',
        normal_sample_id=''
    )

    # Step 3. Write to file
    df_variants = variants_list.to_dataframe()
    df_variants.to_csv(
        get_data_path('hg002_deepvariant.tsv'),
        sep='\t',
        index=False
    )
