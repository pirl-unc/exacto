from .data import get_data_path
from exacto.main import *
from exacto.constants import *


def test_annotate_dna_small_variants_ensembl_gatk4_mutect2():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_gatk4-mutect2.vcf')
    gapped_tsv_file = get_data_path(name='hg38_ucsc_gap_table.txt')
    df_variants = convert_gatk4_mutect2_vcf_to_dataframe(
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
        variant_calling_method=VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2,
        is_tumor_normal_paired=False,
        keep_only_chromosomes=['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'],
        keep_only_filter_values=KEEP_ONLY_FILTER_VALUES,
        min_total_depth=1,
        min_variant_reads_count=1,
        gapped_regions_padding=1
    )

    # Step 3. Annotate
    df_variants_annotated = run_exacto_annotate_genomic_small_variants(
        df_small_variants=df_variants_refined,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95,
        perl_path='',
        annovar_path='',
        annovar_humandb_path='',
        annovar_protocol='',
        annovar_operation='',
        annovar_genome_assembly='',
        annovar_avinput_file='',
        annovar_output_file=''
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
    output_tsv_file = get_data_path('hg002_gatk4-mutect2_annotated_ensembl.tsv')
    df_variants_annotated.to_csv(output_tsv_file, sep='\t', index=False)
