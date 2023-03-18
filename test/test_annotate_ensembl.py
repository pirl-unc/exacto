import pandas as pd
from .data import get_data_path
from exacto.constants import AnnotationSources
from exacto.ensembl import Ensembl
from exacto.variants_list import VariantsList
from exacto.main import run_exacto_annotate


def test_annotate_ensembl():
    # Step 1. Load data
    tsv_file = get_data_path(name='hg002_merged_filtered.tsv')
    variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)

    # Step 2. Load annotation data
    annotation = Ensembl(
        source=AnnotationSources.ENSEMBL,
        release=95,
        species='human'
    )

    # Step 3. Annotate variants
    variants_list = run_exacto_annotate(
        variants_list=variants_list,
        annotation=annotation
    )

    # Step 4. Write to file
    variants_list.to_dataframe().to_csv(
        get_data_path('hg002_merged_filtered_annotated_ensembl.tsv'),
        sep='\t',
        index=False
    )
