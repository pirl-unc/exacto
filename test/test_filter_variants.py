from exacto.constants import VariantFilterQuantifiers, VariantFilterOperators
from exacto.main import run_exacto_filter_variants
from exacto import VariantFilter


def test_filter_cutesv_variants_list(
        cutesv_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=cutesv_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )
    print(variants_list_filtered.size)


def test_filter_pbsv_variants_list1(
        pbsv_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=pbsv_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )
    print(variants_list_filtered.size)


def test_filter_sniffles2_variants_list1(
        sniffles2_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=sniffles2_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )
    print(variants_list_filtered.size)


def test_filter_svim_variants_list1(
        svim_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=svim_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )
    print(variants_list_filtered.size)


def test_filter_deepvariant_variants_list1(
        deepvariant_variants_list,
        hg38_germline_variants_list,
        hg38_excluded_regions_list):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['HG002']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['HG002']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=deepvariant_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )
    print(variants_list_filtered.size)
