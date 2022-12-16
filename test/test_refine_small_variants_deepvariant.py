from .data import get_data_path
from exacto.main import *
from exacto.utilities.vcf_utils import convert_deepvariant_vcf_to_dataframe
from exacto.constants import *
from exacto.default_parameters import *


def test_refine_dna_small_variants_deepvariant():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_deepvariant.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    df_variants = convert_deepvariant_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.PACBIO_HIFI_CCS
    )
    df_gapped_regions = pd.read_csv(gapped_tsv_file, sep='\t')

    # Step 2. Test
    df_variants_refined = run_exacto_refine_genomic_small_variants(
        df_variants=df_variants,
        df_gapped_regions=df_gapped_regions,
        variant_calling_method=VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT,
        tumor_normal_paired=False,
        keep_only_chromosomes=['chr1'],
        keep_only_filter_values=KEEP_ONLY_FILTER_VALUES,
        min_total_coverage=MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
        min_variant_reads_count=MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding=1
    )

    # Step 3. Print output
    print("DataFrame columns:")
    print(df_variants_refined.columns.values.tolist())
    print("%i columns in total" % len(df_variants_refined.columns.values.tolist()))
    print("DataFrame first 5 rows:")
    print(df_variants_refined.head(n=5))
    print("DataFrame first row to dictionary:")
    print(df_variants_refined.iloc[0].to_dict())

