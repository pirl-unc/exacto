from .data import get_data_path
from exacto.main import *
from exacto.constants import *


def test_annotate_dna_structural_variants_ensembl_sniffles2():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_sniffles2.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    germline_sv_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    df_structural_variants = convert_sniffles2_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.PACBIO_HIFI_CCS,
        sample_id='hg002'
    )
    df_gapped_regions = pd.read_csv(gapped_tsv_file, sep='\t')
    df_structural_variants_to_exclude = pd.read_csv(germline_sv_tsv_file, sep='\t')

    # Step 2. Refine
    df_structural_variants_refined = run_exacto_refine_genomic_structural_variants(
        df_structural_variants=df_structural_variants,
        df_structural_variants_to_exclude=df_structural_variants_to_exclude,
        df_gapped_regions=df_gapped_regions,
        variant_calling_method=VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2,
        keep_only_precise_sv=True,
        keep_only_chromosomes=['chr1'],
        keep_only_filter_values=['PASS'],
        min_total_depth=MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count=MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding=GENOME_GAPPED_REGIONS_PADDING,
        exclude_variants_padding=EXCLUDE_SV_PADDING
    )

    # Step 3. Annotate
    df_structural_variants_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_structural_variants_refined,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95
    )

    # Step 4. Print output
    print("First row of DataFrame as dictionary:")
    print(df_structural_variants_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_structural_variants_annotated.columns.values.tolist()))
    print("DataFrame first 5 rows:")
    print(df_structural_variants_annotated.head(n=5))
    print("DataFrame first row to dictionary:")
    print(df_structural_variants_annotated.iloc[0].to_dict())

    # Step 5. Write to file
    output_tsv_file = get_data_path('hg002_sniffles2_annotated_ensembl.tsv')
    df_structural_variants_annotated.to_csv(output_tsv_file, sep='\t', index=False)
