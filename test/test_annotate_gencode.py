import pandas as pd
from .data import get_data_path
from exacto.constants import AnnotationSources
from exacto.gencode import Gencode
from exacto.variants_list import VariantsList
from exacto.main import run_exacto_annotate


def test_annotate_gencode():
    # Step 1. Load data
    tsv_file = get_data_path(name='hg002_merged_filtered.tsv')
    gencode_gtf_file = get_data_path(name='gencode.v41.annotations.gtf')
    variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)

    # Step 2. Load annotation data
    annotation = Gencode(source=AnnotationSources.GENCODE)
    annotation.read_comprehensive_gene_annotation_gtf_file(
        gtf_file=gencode_gtf_file
    )

    # Step 3. Annotate variants
    variants_list = run_exacto_annotate(
        variants_list=variants_list,
        annotation=annotation
    )

    # Step 4. Write to file
    variants_list.to_dataframe().to_csv(
        get_data_path('hg002_merged_filtered_annotated_gencode.tsv'),
        sep='\t',
        index=False
    )
