from .data import get_data_path
from exacto.main import *
from exacto.utilities.vcf_utils import convert_cutesv_vcf_to_dataframe
from exacto.constants import *
from exacto.default_parameters import *


def test_refine_dna_structural_variants_cutesv():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_cutesv.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    germline_sv_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    df_structural_variants = convert_cutesv_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.PACBIO_HIFI_CCS,
        sample_id='hg002'
    )
    df_gapped_regions = pd.read_csv(gapped_tsv_file, sep='\t')
    df_structural_variants_to_exclude = pd.read_csv(germline_sv_tsv_file, sep='\t')

    # Step 2. Test
    df_structural_variants_refined = run_exacto_refine_genomic_structural_variants(
        df_structural_variants=df_structural_variants,
        df_structural_variants_to_exclude=df_structural_variants_to_exclude,
        df_gapped_regions=df_gapped_regions,
        variant_calling_method=VariantCallingMethods.StructuralVariantCallingMethods.CUTESV,
        keep_only_precise_sv=True,
        keep_only_chromosomes=['chr1'],
        keep_only_filter_values=['PASS'],
        min_total_depth=MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count=MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding=GENOME_GAPPED_REGIONS_PADDING,
        exclude_variants_padding=EXCLUDE_SV_PADDING
    )

    # Step 3. Print output
    print("DataFrame columns:")
    print(df_structural_variants_refined.columns.values.tolist())
    print("%i columns in total" % len(df_structural_variants_refined.columns.values.tolist()))
    print("DataFrame first 5 rows:")
    print(df_structural_variants_refined.head(n=5))
    print("DataFrame first row to dictionary:")
    print(df_structural_variants_refined.iloc[0].to_dict())

    # Step 4. Check for errors
    assert len(df_structural_variants_refined.columns.values.tolist()) == \
           len(list(STRUCTURAL_VARIANT_ATTRIBUTES.keys())), \
        "There are %i columns in the refined DataFrame. " \
        "Expected column count is %i." % (
        len(df_structural_variants_refined.columns.values.tolist()),
        len(list(STRUCTURAL_VARIANT_ATTRIBUTES.keys())))

    assert str(df_structural_variants_refined.columns.values.tolist()) == \
           str(list(STRUCTURAL_VARIANT_ATTRIBUTES.keys())), \
        "Expected columns are: %s but found %s in the refined DataFrame." % (
            str(list(STRUCTURAL_VARIANT_ATTRIBUTES.keys())),
            str(df_structural_variants_refined.columns.values.tolist())
        )

    # Step 5. Write to file
    output_tsv_file = get_data_path('hg002_cutesv_refined.tsv')
    df_structural_variants_refined.to_csv(output_tsv_file, sep='\t', index=False)
