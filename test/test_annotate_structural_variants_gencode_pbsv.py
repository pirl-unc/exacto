from .data import get_data_path
from exacto.main import *
from exacto.constants import *
from exacto.variants.vcf import *
from exacto.variants.annotations.gencode import *


def test_annotate_dna_structural_variants_gencode_pbsv():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_pbsv.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    germline_sv_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    gencode_gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
    df_structural_variants = convert_pbsv_vcf_to_dataframe(
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
        variant_calling_method=VariantCallingMethods.StructuralVariantCallingMethods.PBSV,
        keep_only_precise_sv=True,
        keep_only_chromosomes=['chr1'],
        keep_only_filter_values=['PASS'],
        min_total_depth=MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH,
        min_variant_reads_count=MIN_GENOMIC_VARIANT_READS_COUNT,
        gapped_regions_padding=GENOME_GAPPED_REGIONS_PADDING,
        exclude_variants_padding=EXCLUDE_SV_PADDING
    )

    # Step 3. Annotate
    df_gencode_genes, df_gencode_transcripts, df_gencode_exons = read_gencode_gtf_file(
        gencode_gtf_file=gencode_gtf_file
    )
    df_structural_variants_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_structural_variants_refined,
        annotation_source=AnnotationSources.GENCODE,
        df_gencode_genes=df_gencode_genes,
        df_gencode_exons=df_gencode_exons,
        ensembl_release=None
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
    output_tsv_file = get_data_path('hg002_pbsv_annotated_gencode.tsv')
    df_structural_variants_annotated.to_csv(output_tsv_file, sep='\t', index=False)
