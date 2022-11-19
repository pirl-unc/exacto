from .data import get_data_path
from exactolib.main import *
from exactolib.constants import *


def test_annotate_dna_small_variants_ensembl():
    # Step 1. Load data
    sniffles2_tsv_file = get_data_path(name='hg002_sniffles2_refined.tsv')
    pbsv_tsv_file = get_data_path(name='hg002_pbsv_refined.tsv')
    cutesv_tsv_file = get_data_path(name='hg002_cutesv_refined.tsv')
    svim_tsv_file = get_data_path(name='hg002_svim_refined.tsv')
    df_sv_sniffles2 = pd.read_csv(sniffles2_tsv_file, sep='\t')
    df_sv_pbsv = pd.read_csv(pbsv_tsv_file, sep='\t')
    df_sv_cutesv = pd.read_csv(cutesv_tsv_file, sep='\t')
    df_sv_svim = pd.read_csv(svim_tsv_file, sep='\t')

    # Step 2. Test
    df_sv_sniffles2_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_sv_sniffles2,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95
    )
    df_sv_pbsv_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_sv_pbsv,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95
    )
    df_sv_cutesv_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_sv_cutesv,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95
    )
    df_sv_svim_annotated = run_exacto_annotate_genomic_structural_variants(
        df_structural_variants=df_sv_svim,
        annotation_source=AnnotationSources.ENSEMBL,
        df_gencode_genes=None,
        df_gencode_exons=None,
        ensembl_release=95
    )

    # Step 3. Print output
    print("First row of sniffles2 DataFrame as dictionary:")
    print(df_sv_sniffles2_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_sv_sniffles2_annotated.columns.values.tolist()))

    print("First row of pbsv DataFrame as dictionary:")
    print(df_sv_pbsv_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_sv_pbsv_annotated.columns.values.tolist()))

    print("First row of cutesv DataFrame as dictionary:")
    print(df_sv_cutesv_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_sv_cutesv_annotated.columns.values.tolist()))

    print("First row of svim DataFrame as dictionary:")
    print(df_sv_svim_annotated.iloc[0].to_dict())
    print("%i columns in total" % len(df_sv_svim_annotated.columns.values.tolist()))

    # Step 4. Check for errors
    assert 'protein_coding' in df_sv_sniffles2_annotated['ensembl_pos_1_gene_type'].values.tolist(), \
        "The value 'protein_coding' is expected to appear in the column 'ensembl_pos_1_gene_type'"
    assert 'PRKCZ' in df_sv_sniffles2_annotated['ensembl_pos_1_gene_name'].values.tolist(), \
        "The value 'PRKCZ' is expected to appear in the column 'ensembl_pos_1_gene_name'"

