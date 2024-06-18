from .data import get_data_path
from exactolib.main import call_rna_variants


def test_call_rna_variants_1():
    variant_calls = call_rna_variants(
        bam_file=get_data_path(name='bam/hg38_tumor_long_read_rna_1.bam'),
        sample_id='tumor_1',
        min_reads=3,
        min_mapping_quality=20,
        num_threads=1,
        chromosomes=['chr17']
    )
    for variant_call in variant_calls:
        print(variant_call)

