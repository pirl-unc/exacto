from .data import get_data_path
from exactolib.main import *


def test_simulate_dna_variants():
    fasta = pysam.FastaFile(filename=get_data_path(name='hg38_cancer_genes.fa'))
    run_exacto_simulate_variants(
        nucleic_acid_type='dna',
        fasta=fasta
    )
