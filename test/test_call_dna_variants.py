from .data import get_data_path
from exactolib.main import call_dna_variants


def test_call_dna_variants_1():
    variant_calls = call_dna_variants(
        bam_file=get_data_path(name='bam/hg38_tumor_long_read_dna_1.bam'),
        sample_id='tumor_1',
        min_reads=3,
        min_mapping_quality=20,
        num_threads=1,
        chromosomes=['chr7','chr21']
    )
    for variant_call in variant_calls:
        print(variant_call)


def test_call_dna_variants_2():
    variant_calls = call_dna_variants(
        bam_file=get_data_path(name='bam/hg38_tumor_long_read_dna_2.bam'),
        sample_id='tumor_2',
        min_reads=3,
        min_mapping_quality=20,
        num_threads=1,
        chromosomes=['chr17']
    )
    for variant_call in variant_calls:
        print(variant_call)
