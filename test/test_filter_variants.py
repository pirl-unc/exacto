from .data import get_data_path
from exacto.constants import VariantCallingMethods, VariantFilterQuantifiers, VariantFilterOperators
from exacto.default_parameters import MERGE_MAX_NEIGHBOR_DISTANCE
from exacto.genomic_ranges_list import GenomicRangesList
from exacto.main import run_exacto_convert_vcf, run_exacto_filter_variants
from exacto.variant_filter import VariantFilter
from exacto.variants_list import VariantsList


def test_filter_variants():
    # Step 1. Load data
    vcf_file = get_data_path(name='hg002_hg001_gatk4-mutect2.vcf')
    germline_variants_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    excluded_regions_tsv_file = get_data_path(name='hg38_ucsc_gap_table.tsv')
    germline_variants_list = VariantsList.read_tsv_file(tsv_file=germline_variants_tsv_file)
    excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=excluded_regions_tsv_file)

    # Step 2. Convert
    variants_list = run_exacto_convert_vcf(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.GATK4_MUTECT2,
        sequencing_platform='pacbio'
    )

    # Step 3. Filter variants
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['hg002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['hg002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )

    # Step 4. Write to file
    df_variants = variants_list_filtered.to_dataframe()
    df_variants.to_csv(
        get_data_path('hg002_hg001_gatk4-mutect2_filtered.tsv'),
        sep='\t',
        index=False
    )
