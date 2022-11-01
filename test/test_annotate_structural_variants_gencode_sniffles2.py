from .data import get_data_path
from exactolib.main import *
from exactolib.constants import *
from exactolib.utilities.gencode_utils import *


def test_annotate_dna_structural_variants_gencode_sniffles2():
    # Step 1. Load data
    sniffles2_tsv_file = get_data_path(name='hg002_sniffles2_refined.tsv')
    gencode_gtf_file = get_data_path(name='gencode.v41.annotation.gtf')
    df_sv_sniffles2 = pd.read_csv(sniffles2_tsv_file, sep='\t')
    df_gencode_genes, df_gencode_transcripts, df_gencode_exons = read_gencode_gtf_file(
        gencode_gtf_file=gencode_gtf_file
    )

    # Step 2. Test
    df_sv_sniffles2_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_sv_sniffles2,
        annotation_source=AnnotationSources.GENCODE,
        df_gencode_genes=df_gencode_genes,
        df_gencode_exons=df_gencode_exons,
        ensembl_release=None
    )

    # Step 3. Print output
    print("First row of sniffles2 DataFrame as dictionary:")
    print(df_sv_sniffles2_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_sv_sniffles2_annotated.columns.values.tolist()))

    # Step 4. Check for errors
    assert 'lncRNA' in df_sv_sniffles2_annotated['gencode_pos_1_gene_type'].values.tolist(), \
        "The value 'lncRNA' is expected to appear in the column 'gencode_pos_1_gene_type'"
    assert 'ENSG00000230021' in df_sv_sniffles2_annotated['gencode_pos_1_gene_name'].values.tolist(), \
        "The value 'ENSG00000230021' is expected to appear in the column 'gencode_pos_1_gene_name'"

