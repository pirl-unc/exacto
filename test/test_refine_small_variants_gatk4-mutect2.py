from .data import get_data_path
from exacto.main import *
from exacto.variants.vcf import *
from exacto.constants import *
from exacto.default_parameters import *


def test_refine_dna_small_variants_gatk4_mutect2():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_gatk4-mutect2.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    df_variants = convert_gatk4_mutect2_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.ILLUMINA,
        sample_id='hg002',
        tumor_sample_id='hg002'
    )
    df_gapped_regions = pd.read_csv(gapped_tsv_file, sep='\t')

    # Step 2. Test
    df_variants_refined = run_exacto_refine_genomic_small_variants(
        df_variants=df_variants,
        df_gapped_regions=df_gapped_regions,
        df_exclude_snv_indel=None,
        variant_calling_method=VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2,
        is_tumor_normal_paired=False,
        keep_only_chromosomes=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'],
        keep_only_filter_values=KEEP_ONLY_FILTER_VALUES,
        min_total_depth=1,
        min_variant_reads_count=1,
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

    # Step 4. Write to file
    output_tsv_file = get_data_path('hg002_gatk4-mutect2_refined.tsv')
    df_variants_refined.to_csv(output_tsv_file, sep='\t', index=False)
