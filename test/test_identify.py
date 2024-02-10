from .data import get_data_path
from exactolib.main import identify_rna_variants


def test_identify_1():
    variant_calls = identify_rna_variants(
        bam_file=get_data_path(name='hg38_tp53_rna_minimap2_mdtagged_sorted.bam'),
        min_reads=3,
        min_mapping_quality=20,
        num_threads=1,
        chromosomes=['chr17']
    )
    print(len(variant_calls))