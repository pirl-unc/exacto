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
annotating variants using ensembl (pyensembl).
"""


import pyensembl
import pandas as pd
from ..logging import get_logger


logger = get_logger(__name__)


def annotate_variant_using_pyensembl(ensembl: pyensembl.EnsemblRelease,
                                     chromosome: str,
                                     position: int) -> dict:
    """
    Annotates a variant using pyensembl.

    Parameters
    ----------
    ensembl         :   An instance of pyensembl.EnsemblRelease
    chromosome      :   Chromosome.
    position        :   Genomic position.

    Returns
    -------
    Dictionary with the following keys:
    'gene_id'
    'gene_name'
    'gene_type'
    'gene_strand'
    'gene_start'
    'gene_end'
    'region'
    """
    data = {
        'gene_id': '',
        'gene_name': '',
        'gene_type': '',
        'gene_strand': '',
        'gene_start': '',
        'gene_end': '',
        'region': ''
    }

    chromosome = chromosome.replace("chr", "")
    genes = ensembl.genes_at_locus(contig=chromosome, position=position)

    if len(genes) == 0:
        return data
    else:
        gene_ids = []
        gene_names = []
        gene_types = []
        gene_strands = []
        gene_start = []
        gene_end = []

        position_exons = ensembl.exons_at_locus(contig=chromosome, position=position)
        position_exon_ids = []
        for curr_exon in position_exons:
            position_exon_ids.append(curr_exon.exon_id)

        is_exonic = False
        for curr_gene in genes:
            gene_ids.append(curr_gene.gene_id)
            gene_names.append(curr_gene.gene_name)
            gene_types.append(curr_gene.biotype)
            gene_strands.append(curr_gene.strand)
            gene_start.append(str(curr_gene.start))
            gene_end.append(str(curr_gene.end))

            if curr_gene.biotype == 'protein_coding':
                curr_gene_exons = ensembl.exon_ids_of_gene_id(curr_gene.gene_id)
                if len(set(position_exon_ids).intersection(set(curr_gene_exons))) > 0:
                   is_exonic = True

        data['gene_id'] = ','.join(gene_ids)
        data['gene_name'] = ','.join(gene_names)
        data['gene_type'] = ','.join(gene_types)
        data['gene_strand'] = ','.join(gene_strands)
        data['gene_start'] = ','.join(gene_start)
        data['gene_end'] = ','.join(gene_end)

        if is_exonic:
            data['region'] = 'exonic'
        else:
            data['region'] = ''

        return data


def annotate_small_variants_using_pyensembl(df_small_variants: pd.DataFrame,
                                            ensembl_release: int,) -> pd.DataFrame:
    """
    Annotates small variants using ENSEMBL.

    Parameters
    ----------
    df_small_variants       :   DataFrame of structural variants.
                                Expected columns:
                                'chrom'
                                'pos'
    ensembl_release         :   Ensembl release (e.g. 106).

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
    # Step 1. Load Ensembl
    ensembl = pyensembl.EnsemblRelease(ensembl_release)

    # Step 2. Annotate each variant
    data = {
        'id': [],
        'ensembl_pos_region': [],
        'ensembl_pos_gene_id': [],
        'ensembl_pos_gene_name': [],
        'ensembl_pos_gene_type': [],
        'ensembl_pos_gene_strand': [],
        'ensembl_pos_gene_start': [],
        'ensembl_pos_gene_end': []
    }
    for index, row in df_small_variants.iterrows():
        data['id'].append(row['id'])

        # Position annotation
        curr_chrom = row['chrom']
        curr_pos = row['pos']
        curr_annotation = annotate_variant_using_pyensembl(
            ensembl=ensembl,
            chromosome=curr_chrom,
            position=curr_pos
        )
        data['ensembl_pos_region'].append(curr_annotation['region'])
        data['ensembl_pos_gene_id'].append(curr_annotation['gene_id'])
        data['ensembl_pos_gene_name'].append(curr_annotation['gene_name'])
        data['ensembl_pos_gene_type'].append(curr_annotation['gene_type'])
        data['ensembl_pos_gene_strand'].append(curr_annotation['gene_strand'])
        data['ensembl_pos_gene_start'].append(curr_annotation['gene_start'])
        data['ensembl_pos_gene_end'].append(curr_annotation['gene_end'])

    df_annotations = pd.DataFrame(data)
    df_small_variants = pd.merge(df_small_variants, df_annotations, on='id')
    return df_small_variants


def annotate_structural_variants_using_pyensembl(df_structural_variants: pd.DataFrame,
                                                 ensembl_release: int,) -> pd.DataFrame:
    """
    Annotates structural variants using ENSEMBL.

    Parameters
    ----------
    df_structural_variants  :   DataFrame of structural variants.
                                Expected columns:
                                'chr_1'
                                'pos_1'
                                'chr_2'
                                'pos_2'
    ensembl_release         :   Ensembl release (e.g. 106).

    Returns
    -------
    DataFrame with the following columns appended:
    'ensembl_pos_1_region'
    'ensembl_pos_1_gene_id'
    'ensembl_pos_1_gene_name'
    'ensembl_pos_1_gene_type'
    'ensembl_pos_1_gene_strand'
    'ensembl_pos_1_gene_start'
    'ensembl_pos_1_gene_end',
    'ensembl_pos_2_region'
    'ensembl_pos_2_gene_id'
    'ensembl_pos_2_gene_name'
    'ensembl_pos_2_gene_type'
    'ensembl_pos_2_gene_strand'
    'ensembl_pos_2_gene_start'
    'ensembl_pos_2_gene_end'
    """
    # Step 1. Load Ensembl
    ensembl = pyensembl.EnsemblRelease(ensembl_release)

    # Step 2. Annotate each variant
    data = {
        'id': [],
        'ensembl_pos_1_region': [],
        'ensembl_pos_1_gene_id': [],
        'ensembl_pos_1_gene_name': [],
        'ensembl_pos_1_gene_type': [],
        'ensembl_pos_1_gene_strand': [],
        'ensembl_pos_1_gene_start': [],
        'ensembl_pos_1_gene_end': [],
        'ensembl_pos_2_region': [],
        'ensembl_pos_2_gene_id': [],
        'ensembl_pos_2_gene_name': [],
        'ensembl_pos_2_gene_type': [],
        'ensembl_pos_2_gene_strand': [],
        'ensembl_pos_2_gene_start': [],
        'ensembl_pos_2_gene_end': []
    }
    for index, row in df_structural_variants.iterrows():
        data['id'].append(row['id'])

        # Position 1 annotation
        curr_chr_1 = row['chr_1']
        curr_pos_1 = row['pos_1']
        curr_annotation_1 = annotate_variant_using_pyensembl(
            ensembl=ensembl,
            chromosome=curr_chr_1,
            position=curr_pos_1
        )
        data['ensembl_pos_1_region'].append(curr_annotation_1['region'])
        data['ensembl_pos_1_gene_id'].append(curr_annotation_1['gene_id'])
        data['ensembl_pos_1_gene_name'].append(curr_annotation_1['gene_name'])
        data['ensembl_pos_1_gene_type'].append(curr_annotation_1['gene_type'])
        data['ensembl_pos_1_gene_strand'].append(curr_annotation_1['gene_strand'])
        data['ensembl_pos_1_gene_start'].append(curr_annotation_1['gene_start'])
        data['ensembl_pos_1_gene_end'].append(curr_annotation_1['gene_end'])

        # Position 2 annotation
        curr_chr_2 = row['chr_2']
        curr_pos_2 = row['pos_2']
        curr_annotation_2 = annotate_variant_using_pyensembl(
            ensembl=ensembl,
            chromosome=curr_chr_2,
            position=curr_pos_2
        )
        data['ensembl_pos_2_region'].append(curr_annotation_2['region'])
        data['ensembl_pos_2_gene_id'].append(curr_annotation_2['gene_id'])
        data['ensembl_pos_2_gene_name'].append(curr_annotation_2['gene_name'])
        data['ensembl_pos_2_gene_type'].append(curr_annotation_2['gene_type'])
        data['ensembl_pos_2_gene_strand'].append(curr_annotation_2['gene_strand'])
        data['ensembl_pos_2_gene_start'].append(curr_annotation_2['gene_start'])
        data['ensembl_pos_2_gene_end'].append(curr_annotation_2['gene_end'])

    df_annotations = pd.DataFrame(data)
    df_structural_variants = pd.merge(df_structural_variants, df_annotations, on='id')
    return df_structural_variants

