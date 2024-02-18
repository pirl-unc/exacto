from .data import get_data_path
from exactolib.main import call_dna_variants


def test_call_dna_variants_1():
    variant_calls = call_dna_variants(
        bam_file=get_data_path(name='hg38_tp53_tumor_long_read_dna_minimap2_mdtagged_sorted.bam'),
        sample_id='test_001',
        min_reads=3,
        min_mapping_quality=20,
        num_threads=1,
        chromosomes=['chr17']
    )
    for variant_call in variant_calls:
        print(variant_call)


