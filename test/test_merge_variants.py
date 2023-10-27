from .conftest import *
from exacto.default import MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE
from exacto import run_exacto_merge_variants


def test_merge_variants(
        cutesv_variants_list,
        deepvariant_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        svim_variants_list):
    variants_lists = [
        cutesv_variants_list,
        deepvariant_variants_list,
        pbsv_variants_list,
        sniffles2_variants_list,
        svim_variants_list
    ]
    variants_list_merged = run_exacto_merge_variants(
        variants_lists=variants_lists,
        num_threads=1,
        max_neighbor_distance=MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE
    )
    print(variants_list_merged.size)
