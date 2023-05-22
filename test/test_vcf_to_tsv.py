from .conftest import *


def test_cutesv_vcf_to_tsv(cutesv_variants_list):
    print(cutesv_variants_list.variant_ids)


def test_deepvariant_vcf_to_tsv(deepvariant_variants_list):
    print(deepvariant_variants_list.variant_ids)


def test_gatk4_mutect2_vcf_to_tsv(gatk4_mutect2_variants_list):
    print(gatk4_mutect2_variants_list.variant_ids)


def test_pbsv_vcf_to_tsv(pbsv_variants_list):
    print(pbsv_variants_list.variant_ids)


def test_sniffles2_vcf_to_tsv(sniffles2_variants_list):
    print(sniffles2_variants_list.variant_ids)


def test_strelka2_indels_vcf_to_tsv(strelka2_indels_variants_list):
    print(strelka2_indels_variants_list.variant_ids)


def test_strelka2_snvs_vcf_to_tsv(strelka2_snvs_variants_list):
    print(strelka2_snvs_variants_list.variant_ids)


def test_svim_vcf_to_tsv(svim_variants_list):
    print(svim_variants_list.variant_ids)

