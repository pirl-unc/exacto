from exactolib.main import identify_rna_variants


def test_identify_1():
    variant_calls = identify_rna_variants(
        bam_file='/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/hg38_tp53_rna_minimap2_mdtagged_sorted.bam',
        min_reads=1,
        min_mapping_quality=20,
        num_threads=2,
        chromosomes=['chr17']
    )
    print(len(variant_calls))