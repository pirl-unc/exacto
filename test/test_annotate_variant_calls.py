from .conftest import *
from .data import get_data_path
from exacto.ensembl import Ensembl
from exacto.gencode import Gencode
from exacto.main import run_exacto_annotate_variants_list


"""Structural Variants"""
def test_annotate_sniffles2_variants_list_by_ensembl(
        sniffles2_variants_list,
):
    ensembl = Ensembl(release=95, species='human')
    variants_list_annotated = run_exacto_annotate_variants_list(
        variants_list=sniffles2_variants_list,
        annotation_db=ensembl
    )
    print(variants_list_annotated.variant_ids)

# def test_annotate_sniffles2_variants_list_by_gencode(
#         sniffles2_variants_list,
# ):
#     gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
#     genome_fasta_file = get_data_path(name='hg38.fa')
#     gencode = Gencode(
#         gtf_file=gtf_file,
#         genome_fasta_file=genome_fasta_file,
#         version='v41',
#         genome='hg38',
#     )
#     variants_list_annotated = run_exacto_annotate_variants_list(
#         variants_list=sniffles2_variants_list,
#         annotation_db=gencode
#     )
#     print(variants_list_annotated.variant_ids)
