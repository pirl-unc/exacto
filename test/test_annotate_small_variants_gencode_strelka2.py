from .data import get_data_path
from exacto.main import *
from exacto.variants.vcf import *
from exacto.variants.annotations.gencode import *
from exacto.constants import *


def test_annotate_dna_small_variants_gencode_strelka2():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_strelka2.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    gencode_gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
    df_variants = convert_strelka2_germline_vcf_to_dataframe(
        vcf_file=vcf_file,
        sequencing_platform=SequencingPlatforms.PACBIO_HIFI_CCS,
        sample_id='hg002',
        tumor_sample_id='hg002'
    )
    df_gapped_regions = pd.read_csv(gapped_tsv_file, sep='\t')

    # Step 2. Refine
    df_variants_refined = run_exacto_refine_genomic_small_variants(
        df_variants=df_variants,
        df_gapped_regions=df_gapped_regions,
        df_exclude_snv_indel=None,
        variant_calling_method=VariantCallingMethods.SmallVariantCallingMethods.STRELKA2_GERMLINE,
        is_tumor_normal_paired=False,
        keep_only_chromosomes=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'],
        keep_only_filter_values=KEEP_ONLY_FILTER_VALUES,
        min_total_depth=1,
        min_variant_reads_count=1,
        gapped_regions_padding=1
    )

    # Step 3. Annotate
    df_gencode_genes, df_gencode_transcripts, df_gencode_exons = read_gencode_gtf_file(
        gencode_gtf_file=gencode_gtf_file
    )
    df_variants_annotated = run_exacto_annotate_genomic_small_variants(
        df_small_variants=df_variants_refined,
        annotation_source=AnnotationSources.GENCODE,
        df_gencode_genes=df_gencode_genes,
        df_gencode_exons=df_gencode_exons,
    )

    # Step 4. Print output
    print("First row of DataFrame as dictionary:")
    print(df_variants_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_variants_annotated.columns.values.tolist()))
    print("DataFrame first 5 rows:")
    print(df_variants_annotated.head(n=5))
    print("DataFrame first row to dictionary:")
    print(df_variants_annotated.iloc[0].to_dict())

    # Step 5. Write to file
    output_tsv_file = get_data_path('hg002_strelka2_annotated_gencode.tsv')
    df_variants_annotated.to_csv(output_tsv_file, sep='\t', index=False)
