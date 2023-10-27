from .conftest import *


def test_cutesv_vcf2tsv(cutesv_variants_list):
    print(cutesv_variants_list.size)


def test_deepvariant_vcf2tsv(deepvariant_variants_list):
    print(deepvariant_variants_list.size)


def test_gatk4_mutect2_vcf2tsv(gatk4_mutect2_variants_list):
    print(gatk4_mutect2_variants_list.size)


def test_pbsv_vcf2tsv(pbsv_variants_list):
    print(pbsv_variants_list.size)


def test_sniffles2_vcf2tsv(sniffles2_variants_list):
    print(sniffles2_variants_list.size)


def test_strelka2_indels_vcf2tsv(strelka2_indels_variants_list):
    print(strelka2_indels_variants_list.size)


def test_strelka2_snvs_vcf2tsv(strelka2_snvs_variants_list):
    print(strelka2_snvs_variants_list.size)


def test_svim_vcf2tsv(svim_variants_list):
    print(svim_variants_list.size)

