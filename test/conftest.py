import pytest
from .data import get_data_path
from exacto.constants import *
from exacto.genomic_ranges_list import GenomicRangesList
from exacto.main import run_exacto_vcf_to_tsv
from exacto.variants_list import VariantsList


"""Structural Variants Lists"""
@pytest.fixture
def cutesv_variants_list():
    vcf_file = get_data_path(name='hg002_cutesv.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.CUTESV,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def pbsv_variants_list():
    vcf_file = get_data_path(name='hg002_pbsv.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.PBSV,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def sniffles2_variants_list():
    vcf_file = get_data_path(name='hg002_sniffles2.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SNIFFLES2,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def svim_variants_list():
    vcf_file = get_data_path(name='hg002_svim.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.SVIM,
        sequencing_platform='pacbio'
    )
    return variants_list


"""Small Variants Lists"""
@pytest.fixture
def deepvariant_variants_list():
    vcf_file = get_data_path(name='hg002_deepvariant.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.DEEPVARIANT,
        sequencing_platform='pacbio'
    )
    return variants_list

@pytest.fixture
def gatk4_mutect2_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_gatk4-mutect2.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.GATK4_MUTECT2,
        sequencing_platform='illumina'
    )
    return variants_list

@pytest.fixture
def strelka2_indels_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_strelka2_indels.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.STRELKA2,
        sequencing_platform='illumina'
    )
    return variants_list

@pytest.fixture
def strelka2_snvs_variants_list():
    vcf_file = get_data_path(name='hg002_hg001_strelka2_snvs.vcf')
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=vcf_file,
        source_id='hg002',
        variant_calling_method=VariantCallingMethods.STRELKA2,
        sequencing_platform='illumina'
    )
    return variants_list


"""Germline Variants List"""
@pytest.fixture
def germline_variants_list():
    germline_variants_tsv_file = get_data_path(name='audano_et_al_cell_2019_sv_list.tsv')
    germline_variants_list = VariantsList.read_tsv_file(tsv_file=germline_variants_tsv_file)
    return germline_variants_list


"""Excluded Regions List"""
@pytest.fixture
def excluded_regions_list():
    excluded_regions_tsv_file = get_data_path(name='hg38_ucsc_gap_table.tsv')
    excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=excluded_regions_tsv_file)
    return excluded_regions_list

