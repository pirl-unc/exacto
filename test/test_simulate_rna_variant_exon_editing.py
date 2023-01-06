from .data import get_data_path
from exacto.simulation.rna.match import Match
from exacto.simulation.rna.deletion import Deletion
from exacto.simulation.rna.insertion import Insertion
from exacto.simulation.rna.substitution import Substitution
from exacto.constants import *


def test_simulate_rna_variant_exon_editing():
    ref_exon = Match(
        gene_id='gene001',
        transcript_id='transcript001',
        exon_id='exon001',
        exon_number=1,
        strand='+',
        chrom='chr1',
        start=651,
        end=660,
        length=10,
        sequence='ATCGCCATTC',
    )
    print(ref_exon)

    del_exon = Deletion(
        ref_exon=ref_exon,
        del_start=652,
        del_end=653
    )
    print(del_exon)

    ins_exon = Insertion(
        ref_exon=del_exon,
        ins_pos=655,
        ins_sequence='TT'
    )
    print(ins_exon)

    snv_exon = Substitution(
        ref_exon=ins_exon,
        snv_pos=657,
        snv_alt='G'
    )
    print(snv_exon)

    del_exon = Deletion(
        ref_exon=snv_exon,
        del_start=658,
        del_end=659
    )
    print(del_exon)
