#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions related to
running pyensembl.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import pyensembl
import pandas as pd
import subprocess as sp
import os
from exactolib.logging import get_logger


logger = get_logger(__name__)


def annotate_variant(ensembl, chromosome, position):
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
        return __
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


def annotate_using_pyensembl(tsv_file,
                             ensembl_release,
                             output_tsv_file):
    """
    Annotates a TSV file using pyensembl

    Returns
    -------

    """
    # Step 1. Load Ensembl
    ensembl = pyensembl.EnsemblRelease(ensembl_release)

    # Step 2. Read TSV file
    df_variants = pd.read_csv(tsv_file, sep='\t')

    # Step 3. Annotate each variant
    data = {
        'id': [],
        'pos_1_region': [],
        'pos_1_gene_id': [],
        'pos_1_gene_name': [],
        'pos_1_gene_type': [],
        'pos_1_gene_strand': [],
        'pos_1_gene_start': [],
        'pos_1_gene_end': [],
        'pos_2_region': [],
        'pos_2_gene_id': [],
        'pos_2_gene_name': [],
        'pos_2_gene_type': [],
        'pos_2_gene_strand': [],
        'pos_2_gene_start': [],
        'pos_2_gene_end': []
    }
    for index, row in df_variants.iterrows():
        data['id'].append(row['id'])

        # Position 1 annotation
        curr_chr_1 = row['chr_1']
        curr_pos_1 = row['pos_1']
        curr_annotation_1 = annotate_variant(ensembl=ensembl,
                                             chromosome=curr_chr_1,
                                             position=curr_pos_1)
        data['pos_1_region'].append(curr_annotation_1['region'])
        data['pos_1_gene_id'].append(curr_annotation_1['gene_id'])
        data['pos_1_gene_name'].append(curr_annotation_1['gene_name'])
        data['pos_1_gene_type'].append(curr_annotation_1['gene_type'])
        data['pos_1_gene_strand'].append(curr_annotation_1['gene_strand'])
        data['pos_1_gene_start'].append(curr_annotation_1['gene_start'])
        data['pos_1_gene_end'].append(curr_annotation_1['gene_end'])

        # Position 2 annotation
        curr_chr_2 = row['chr_2']
        curr_pos_2 = row['pos_2']
        curr_annotation_2 = annotate_variant(ensembl=ensembl,
                                             chromosome=curr_chr_2,
                                             position=curr_pos_2)
        data['pos_2_region'].append(curr_annotation_2['region'])
        data['pos_2_gene_id'].append(curr_annotation_2['gene_id'])
        data['pos_2_gene_name'].append(curr_annotation_2['gene_name'])
        data['pos_2_gene_type'].append(curr_annotation_2['gene_type'])
        data['pos_2_gene_strand'].append(curr_annotation_2['gene_strand'])
        data['pos_2_gene_start'].append(curr_annotation_2['gene_start'])
        data['pos_2_gene_end'].append(curr_annotation_2['gene_end'])

    df_annotations = pd.DataFrame(data)
    df = pd.merge(df_variants, df_annotations, on='id')
    df.to_csv(output_tsv_file, sep='\t', index=False)

    return 0

