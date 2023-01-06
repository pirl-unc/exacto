from .data import get_data_path
from exacto.main import *
from exacto.constants import *


def test_convert_deepvariant_vcf():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_sniffles2.vcf')

    # Step 2. Convert
    df_variants = convert_sniffles2_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.PACBIO_HIFI_CCS,
        sample_id='hg002'
    )

    # Step 3. Write to file
    output_tsv_file = get_data_path('hg002_sniffles2.tsv')
    df_variants.to_csv(output_tsv_file, sep='\t', index=False)
