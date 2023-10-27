import pandas as pd
from exacto import VariantsList
from exacto import run_exacto_vcf2tsv
from exacto import run_exacto_merge_variants
from exacto.constants import VariantCallingMethods
from exacto.default import MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE
from exacto.constants import VariantFilterQuantifiers, VariantFilterOperators
from exacto.main import run_exacto_filter_variants
from exacto import VariantFilter, GenomicRangesList


VCF_FILE_1 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg002_sniffles2.vcf"
VCF_FILE_2 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg002_cutesv.vcf"
VCF_FILE_3 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg002_svim.vcf"
VCF_FILE_4 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg002_deepvariant.vcf"
VCF_FILE_5 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg002_pbsv.vcf"

TSV_FILE_1 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/audano_et_al_cell_2019_sv_list.tsv"
TSV_FILE_2 = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg38_ucsc_gap_table.tsv"


if __name__ == "__main__":
    variants_list_1 = run_exacto_vcf2tsv(
        vcf_file=VCF_FILE_1,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SNIFFLES2,
        sequencing_platform='pacbio'
    )
    variants_list_2 = run_exacto_vcf2tsv(
        vcf_file=VCF_FILE_2,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.CUTESV,
        sequencing_platform='pacbio'
    )
    variants_list_3 = run_exacto_vcf2tsv(
        vcf_file=VCF_FILE_3,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SVIM,
        sequencing_platform='pacbio'
    )
    variants_list_4 = run_exacto_vcf2tsv(
        vcf_file=VCF_FILE_4,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.DEEPVARIANT,
        sequencing_platform='pacbio'
    )
    variants_list_5 = run_exacto_vcf2tsv(
        vcf_file=VCF_FILE_5,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.PBSV,
        sequencing_platform='pacbio'
    )

    hg38_germline_variants_list = VariantsList.read_tsv_file(tsv_file=TSV_FILE_1)
    hg38_excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=TSV_FILE_2)

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
        variants_list=variants_list_5,
        variant_filters=variant_filters,
        excluded_variants_list=hg38_germline_variants_list,
        excluded_regions_list=hg38_excluded_regions_list,
        num_threads=1
    )

    # variants_lists = [
    #     variants_list_1,
    #     variants_list_2,
    #     variants_list_3,
    #     variants_list_4,
    #     variants_list_5
    # ]
    # variants_list_merged = run_exacto_merge_variants(
    #     variants_lists=variants_lists,
    #     num_threads=1,
    #     max_neighbor_distance=MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE
    # )
    # print(variants_list_merged.size)
