from .conftest import *
from exacto.default_parameters import MERGE_MAX_NEIGHBOR_DISTANCE
from exacto.main import run_exacto_merge_variant_calls


def test_merge_structural_variant_calls(
        cutesv_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        svim_variants_list
):
    variants_lists = [
        cutesv_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        svim_variants_list
    ]
    variants_list_merged = run_exacto_merge_variant_calls(
        variants_lists=variants_lists,
        max_neighbor_distance=MERGE_MAX_NEIGHBOR_DISTANCE
    )
    print(variants_list_merged.variant_ids)


def test_merge_small_variant_calls(
        gatk4_mutect2_variants_list,
        strelka2_snvs_variants_list,
        strelka2_indels_variants_list
):
    variants_lists = [
        gatk4_mutect2_variants_list,
        strelka2_snvs_variants_list,
        strelka2_indels_variants_list
    ]
    variants_list_merged = run_exacto_merge_variant_calls(
        variants_lists=variants_lists,
        max_neighbor_distance=MERGE_MAX_NEIGHBOR_DISTANCE
    )
    print(variants_list_merged.variant_ids)


def test_merge_all_variant_calls(
        cutesv_variants_list,
        deepvariant_variants_list,
        gatk4_mutect2_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        strelka2_snvs_variants_list,
        strelka2_indels_variants_list,
        svim_variants_list
):
    variants_lists = [
        cutesv_variants_list,
        deepvariant_variants_list,
        gatk4_mutect2_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        strelka2_snvs_variants_list,
        strelka2_indels_variants_list,
        svim_variants_list
    ]
    variants_list_merged = run_exacto_merge_variant_calls(
        variants_lists=variants_lists,
        max_neighbor_distance=MERGE_MAX_NEIGHBOR_DISTANCE
    )
    print(variants_list_merged.variant_ids)
