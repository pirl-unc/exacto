from .conftest import *
from exacto.constants import VariantFilterQuantifiers, VariantFilterOperators
from exacto.main import run_exacto_filter_variants
from exacto.variant_filter import VariantFilter


"""Structural Variants"""
def test_filter_cutesv_variants_list(
        cutesv_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=cutesv_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_pbsv_variants_list(
        pbsv_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=pbsv_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_sniffles2_variants_list(
        sniffles2_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=sniffles2_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_svim_variants_list(
        svim_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=svim_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)


"""Small Variants"""
def test_filter_deepvariant_variants_list(
        deepvariant_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=deepvariant_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_gatk4_mutect2_variants_list(
        gatk4_mutect2_variants_list,
        germline_variants_list,
        excluded_regions_list
):
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
        variants_list=gatk4_mutect2_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_strelka2_indels_variants_list(
        strelka2_indels_variants_list,
        germline_variants_list,
        excluded_regions_list
):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['TUMOR']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['TUMOR']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=strelka2_indels_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

def test_filter_strelka2_snvs_variants_list(
        strelka2_snvs_variants_list,
        germline_variants_list,
        excluded_regions_list
):
    variant_filters = []
    variant_filter_1 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='alternate_allele_read_count',
        operator=VariantFilterOperators.GREATER_THAN_OR_EQUAL_TO,
        value=2,
        sample_ids=['TUMOR']
    )
    variant_filter_2 = VariantFilter(
        quantifier=VariantFilterQuantifiers.ALL,
        attribute='filter',
        operator=VariantFilterOperators.EQUALS,
        value='PASS',
        sample_ids=['TUMOR']
    )
    variant_filters.append(variant_filter_1)
    variant_filters.append(variant_filter_2)
    variants_list_filtered = run_exacto_filter_variants(
        variants_list=strelka2_snvs_variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=germline_variants_list,
        excluded_regions_list=excluded_regions_list
    )
    print(variants_list_filtered.variant_ids)

