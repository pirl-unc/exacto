# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
The purpose of this python3 script is to implement functions related to
annotating variants using GENCODE.
"""


import pandas as pd
from ..logging import get_logger


logger = get_logger(__name__)


def annotate_small_variants_using_gencode(df_small_variants: pd.DataFrame,
                                          df_gencode_genes: pd.DataFrame,
                                          df_gencode_exons: pd.DataFrame) -> pd.DataFrame:
    """
    Annotates small variants using GENCODE.

    Parameters
    ----------
    df_small_variants       :   DataFrame of small variants.
                                Expected columns:
                                'id'
                                'chrom'
                                'pos'
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'

    Returns
    -------
    DataFrame with the following columns appended:
    'ensembl_pos_region'
    'ensembl_pos_gene_id'
    'ensembl_pos_gene_name'
    'ensembl_pos_gene_type'
    'ensembl_pos_gene_strand'
    'ensembl_pos_gene_start'
    'ensembl_pos_gene_end'
    """
    # Step 1. Annotate each variant
    data = {
        'id': [],
        'gencode_pos_region': [],
        'gencode_pos_gene_id': [],
        'gencode_pos_gene_name': [],
        'gencode_pos_gene_type': [],
        'gencode_pos_gene_strand': [],
        'gencode_pos_gene_start': [],
        'gencode_pos_gene_end': [],
        'gencode_pos_exon_id': [],
        'gencode_pos_exon_number': []
    }
    for index, row in df_small_variants.iterrows():
        data['id'].append(row['id'])

        # Position 1 annotation
        curr_chrom = row['chrom']
        curr_pos = row['pos']
        curr_pos_region = ''
        curr_pos_gene_id = ''
        curr_pos_gene_name = ''
        curr_pos_gene_type = ''
        curr_pos_gene_strand = ''
        curr_pos_gene_start = ''
        curr_pos_gene_end = ''
        curr_pos_exon_id = ''
        curr_pos_exon_number = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chrom) &
            (df_gencode_genes['gene_start'] <= curr_pos) &
            (df_gencode_genes['gene_end'] >= curr_pos),:
        ]
        if len(df_gencode_genes_matched) > 0:
             df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chrom) &
                (df_gencode_exons['exon_start'] <= curr_pos) &
                (df_gencode_exons['exon_end'] >= curr_pos), :
             ]
             if len(df_gencode_exons_matched) > 0:
                 curr_pos_region = 'exonic'
                 curr_pos_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                 curr_pos_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
             else:
                 curr_pos_region = 'intronic'
             curr_pos_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
             curr_pos_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
             curr_pos_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
             curr_pos_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
             curr_pos_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
             curr_pos_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))

        else:
            curr_pos_region = 'intergenic'
        data['gencode_pos_region'].append(curr_pos_region)
        data['gencode_pos_gene_id'].append(curr_pos_gene_id)
        data['gencode_pos_gene_name'].append(curr_pos_gene_name)
        data['gencode_pos_gene_type'].append(curr_pos_gene_type)
        data['gencode_pos_gene_strand'].append(curr_pos_gene_strand)
        data['gencode_pos_gene_start'].append(curr_pos_gene_start)
        data['gencode_pos_gene_end'].append(curr_pos_gene_end)
        data['gencode_pos_exon_id'].append(curr_pos_exon_id)
        data['gencode_pos_exon_number'].append(curr_pos_exon_number)

    df_annotations = pd.DataFrame(data)
    df_small_variants = pd.merge(df_small_variants, df_annotations, on='id')
    return df_small_variants


def annotate_structural_variants_using_gencode(df_structural_variants: pd.DataFrame,
                                               df_gencode_genes: pd.DataFrame,
                                               df_gencode_exons: pd.DataFrame) -> pd.DataFrame:
    """
    Annotates structural variants using GENCODE.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of structural variants.
                                Expected columns:
                                'id'
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
                                'sv_type' (DEL, INS, INV, DUP, BND or TRA)
    df_gencode_genes        :   DataFrame of GENCODE genes.
                                Expected columns:
                                'gene_id'
                                'gene_name'
                                'gene_type'
                                'gene_chrom'
                                'gene_start'
                                'gene_end'
                                'gene_strand'
                                'level'
                                'transcripts_count'
    df_gencode_exons        :   DataFrame of GENCODE exons.
                                Expected columns:
                                'gene_id'
                                'transcript_id'
                                'exon_id'
                                'exon_number'
                                'exon_chrom'
                                'exon_start'
                                'exon_end'

    Returns
    -------
    DataFrame with the following columns appended:
    'ensembl_pos_1_region'
    'ensembl_pos_1_gene_id'
    'ensembl_pos_1_gene_name'
    'ensembl_pos_1_gene_type'
    'ensembl_pos_1_gene_strand'
    'ensembl_pos_1_gene_start'
    'ensembl_pos_1_gene_end'
    'ensembl_pos_2_region'
    'ensembl_pos_2_gene_id'
    'ensembl_pos_2_gene_name'
    'ensembl_pos_2_gene_type'
    'ensembl_pos_2_gene_strand'
    'ensembl_pos_2_gene_start'
    'ensembl_pos_2_gene_end'
    """
    # Step 1. Annotate each variant
    data = {
        'id': [],
        'gencode_pos_1_region': [],
        'gencode_pos_1_gene_id': [],
        'gencode_pos_1_gene_name': [],
        'gencode_pos_1_gene_type': [],
        'gencode_pos_1_gene_strand': [],
        'gencode_pos_1_gene_start': [],
        'gencode_pos_1_gene_end': [],
        'gencode_pos_1_exon_id': [],
        'gencode_pos_1_exon_number': [],
        'gencode_pos_2_region': [],
        'gencode_pos_2_gene_id': [],
        'gencode_pos_2_gene_name': [],
        'gencode_pos_2_gene_type': [],
        'gencode_pos_2_gene_strand': [],
        'gencode_pos_2_gene_start': [],
        'gencode_pos_2_gene_end': [],
        'gencode_pos_2_exon_id': [],
        'gencode_pos_2_exon_number': []
    }
    for index, row in df_structural_variants.iterrows():
        data['id'].append(row['id'])

        # Position 1 annotation
        curr_chr_1 = row['chr_1']
        curr_pos_1 = row['pos_1']
        curr_pos_1_region = ''
        curr_pos_1_exon_id = ''
        curr_pos_1_exon_number = ''
        curr_pos_1_gene_id = ''
        curr_pos_1_gene_name = ''
        curr_pos_1_gene_type = ''
        curr_pos_1_gene_strand = ''
        curr_pos_1_gene_start = ''
        curr_pos_1_gene_end = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chr_1) &
            (df_gencode_genes['gene_start'] <= curr_pos_1) &
            (df_gencode_genes['gene_end'] >= curr_pos_1),:
        ]
        if len(df_gencode_genes_matched) > 0:
             df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chr_1) &
                (df_gencode_exons['exon_start'] <= curr_pos_1) &
                (df_gencode_exons['exon_end'] >= curr_pos_1), :
            ]
             if len(df_gencode_exons_matched) > 0:
                 curr_pos_1_region = 'exonic'
                 curr_pos_1_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                 curr_pos_1_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
             else:
                 curr_pos_1_region = 'intronic'
             curr_pos_1_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
             curr_pos_1_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
             curr_pos_1_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
             curr_pos_1_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
             curr_pos_1_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
             curr_pos_1_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))

        else:
            curr_pos_1_region = 'intergenic'
        data['gencode_pos_1_region'].append(curr_pos_1_region)
        data['gencode_pos_1_gene_id'].append(curr_pos_1_gene_id)
        data['gencode_pos_1_gene_name'].append(curr_pos_1_gene_name)
        data['gencode_pos_1_gene_type'].append(curr_pos_1_gene_type)
        data['gencode_pos_1_gene_strand'].append(curr_pos_1_gene_strand)
        data['gencode_pos_1_gene_start'].append(curr_pos_1_gene_start)
        data['gencode_pos_1_gene_end'].append(curr_pos_1_gene_end)
        data['gencode_pos_1_exon_id'].append(curr_pos_1_exon_id)
        data['gencode_pos_1_exon_number'].append(curr_pos_1_exon_number)

        # Position 2 annotation
        curr_chr_2 = row['chr_2']
        curr_pos_2 = row['pos_2']
        curr_pos_2_region = ''
        curr_pos_2_exon_id = ''
        curr_pos_2_exon_number = ''
        curr_pos_2_gene_id = ''
        curr_pos_2_gene_name = ''
        curr_pos_2_gene_type = ''
        curr_pos_2_gene_strand = ''
        curr_pos_2_gene_start = ''
        curr_pos_2_gene_end = ''
        df_gencode_genes_matched = df_gencode_genes.loc[
            (df_gencode_genes['gene_chrom'] == curr_chr_2) &
            (df_gencode_genes['gene_start'] <= curr_pos_2) &
            (df_gencode_genes['gene_end'] >= curr_pos_2), :
        ]
        if len(df_gencode_genes_matched) > 0:
            df_gencode_exons_matched = df_gencode_exons.loc[
                (df_gencode_exons['exon_chrom'] == curr_chr_2) &
                (df_gencode_exons['exon_start'] <= curr_pos_2) &
                (df_gencode_exons['exon_end'] >= curr_pos_2), :
            ]
            if len(df_gencode_exons_matched) > 0:
                curr_pos_2_region = 'exonic'
                curr_pos_2_exon_id = ','.join(df_gencode_exons_matched['exon_id'].values.tolist())
                curr_pos_2_exon_number = ','.join(map(str, df_gencode_exons_matched['exon_number'].values.tolist()))
            else:
                curr_pos_2_region = 'intronic'
            curr_pos_2_gene_id = ','.join(df_gencode_genes_matched['gene_id'].values.tolist())
            curr_pos_2_gene_name = ','.join(df_gencode_genes_matched['gene_name'].values.tolist())
            curr_pos_2_gene_type = ','.join(df_gencode_genes_matched['gene_type'].values.tolist())
            curr_pos_2_gene_strand = ','.join(df_gencode_genes_matched['gene_strand'].values.tolist())
            curr_pos_2_gene_start = ','.join(map(str, df_gencode_genes_matched['gene_start'].values.tolist()))
            curr_pos_2_gene_end = ','.join(map(str, df_gencode_genes_matched['gene_end'].values.tolist()))
        else:
            curr_pos_2_region = 'intergenic'
        data['gencode_pos_2_region'].append(curr_pos_2_region)
        data['gencode_pos_2_gene_id'].append(curr_pos_2_gene_id)
        data['gencode_pos_2_gene_name'].append(curr_pos_2_gene_name)
        data['gencode_pos_2_gene_type'].append(curr_pos_2_gene_type)
        data['gencode_pos_2_gene_strand'].append(curr_pos_2_gene_strand)
        data['gencode_pos_2_gene_start'].append(curr_pos_2_gene_start)
        data['gencode_pos_2_gene_end'].append(curr_pos_2_gene_end)
        data['gencode_pos_2_exon_id'].append(curr_pos_2_exon_id)
        data['gencode_pos_2_exon_number'].append(curr_pos_2_exon_number)

    df_annotations = pd.DataFrame(data)
    df_structural_variants = pd.merge(df_structural_variants, df_annotations, on='id')
    return df_structural_variants

