#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions
related to merging SV TSV files.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import pandas as pd
from exactolib.logging import get_logger


logger = get_logger(__name__)


def make_sorter(sorted_list):
    sort_order = { k:v for k, v in zip(sorted_list, range(len(sorted_list))) }
    return lambda s: s.map(lambda x: sort_order[x])


def merge_sv_tsv_files(tsv_files,
                       methods_priority_list,
                       output_merged_tsv_file,
                       output_merged_deduped_tsv_file,
                       max_cluster_distance=10):
    """
    Merges structural variant TSV files.

    Parameters
    ----------
    tsv_files                   : List of TSV files including full paths.
                                  The expected column headers of each TSV file are:
                                  chr_1,
                                  pos_1,
                                  chr_2,
                                  pos_2,
                                  sv_type,
                                  sv_size,
                                  method,
                                  insertion_sequence,
                                  strand,
                                  reference_reads_count,
                                  variant_reads_count,
                                  total_coverage,
                                  variant_allele_fraction,
                                  genotype,
                                  pos_1_region,
                                  pos_1_gene_id,
                                  pos_1_gene_name,
                                  pos_1_gene_type,
                                  pos_1_gene_strand,
                                  pos_1_gene_start,
                                  pos_1_gene_end,
                                  pos_2_region,
                                  pos_2_gene_id,
                                  pos_2_gene_name,
                                  pos_2_gene_type,
                                  pos_2_gene_strand,
                                  pos_2_gene_start,
                                  pos_2_gene_end,
                                  read_ids (optional; expected from Sniffles VCF)
    methods_priority_list       : List of method used to select a SV when there are duplicates.
    output_tsv_file             : Full path of output TSV file.
    max_cluster_distance        : Maximum distance between methods to cluster (default: 10).

    Returns
    -------
    dataframe
    """
    expected_columns = [
        'chr_1',
        'pos_1',
        'chr_2',
        'pos_2',
        'sv_type',
        'sv_size',
        'method',
        'insertion_sequence',
        'strand',
        'reference_reads_count',
        'variant_reads_count',
        'total_coverage',
        'variant_allele_fraction',
        'genotype',
        'read_ids',
        'pos_1_region',
        'pos_1_gene_id',
        'pos_1_gene_name',
        'pos_1_gene_type',
        'pos_1_gene_strand',
        'pos_1_gene_start',
        'pos_1_gene_end',
        'pos_2_region',
        'pos_2_gene_id',
        'pos_2_gene_name',
        'pos_2_gene_type',
        'pos_2_gene_strand',
        'pos_2_gene_start',
        'pos_2_gene_end'
    ]
    df_all = pd.DataFrame()
    logger.info('Started reading the TSV files')
    i = 0
    for curr_tsv_file in tsv_files:
        df_curr = pd.read_csv(curr_tsv_file,
                              sep='\t',
                              low_memory=False,
                              memory_map=True)
        if 'insertion_sequence' not in df_curr.columns.values.tolist():
            df_curr['insertion_sequence'] = '' * len(df_curr)
        if 'strand' not in df_curr.columns.values.tolist():
            df_curr['strand'] = '' * len(df_curr)
        if 'read_ids' not in df_curr.columns.values.tolist():
            df_curr['read_ids'] = '' * len(df_curr)
        df_curr = df_curr.loc[:, expected_columns]
        df_all = df_all.append(df_curr)
        i += 1
    df_all['record_id'] = range(0, len(df_all))
    df_all.to_csv(output_merged_tsv_file, sep='\t', index=False)

    logger.info('Finished reading the TSV files')

    # Merge SV calls
    df_merged = pd.DataFrame()
    recorded_ids = set()
    n = len(df_all)
    logger.info("%i SV calls to iterate" % n)
    count = 0
    for index, row in df_all.iterrows():
        if row['record_id'] in recorded_ids:
            continue
        curr_chr_1 = str(row['chr_1'])
        curr_pos_1 = int(row['pos_1'])
        curr_chr_2 = str(row['chr_2'])
        curr_pos_2 = int(row['pos_2'])
        curr_sv_type = str(row['sv_type'])

        df_matched = df_all.loc[
            (
                (df_all.sv_type == curr_sv_type) &
                ((df_all.chr_1 == curr_chr_1) & (df_all.chr_2 == curr_chr_2) &
                 (df_all.pos_1 >= curr_pos_1 - max_cluster_distance) &
                 (df_all.pos_1 <= curr_pos_1 + max_cluster_distance) &
                 (df_all.pos_2 >= curr_pos_2 - max_cluster_distance) &
                 (df_all.pos_2 <= curr_pos_2 + max_cluster_distance)) |
                ((df_all.chr_1 == curr_chr_2) & (df_all.chr_2 == curr_chr_1) &
                 (df_all.pos_1 >= curr_pos_2 - max_cluster_distance) &
                 (df_all.pos_1 <= curr_pos_2 + max_cluster_distance) &
                 (df_all.pos_2 >= curr_pos_1 - max_cluster_distance) &
                 (df_all.pos_2 <= curr_pos_1 + max_cluster_distance))
            ),:]

        df_matched = df_matched.sort_values('method', key=make_sorter(methods_priority_list))

        for curr_record_id in df_matched.record_id.values.tolist():
            recorded_ids.add(curr_record_id)

        # Record data
        methods = df_matched.method.unique()
        methods_count = len(methods)
        df_curr = df_matched.iloc[0].copy()
        df_curr['methods'] = ','.join(methods)
        df_curr['methods_count'] = methods_count
        df_merged = df_merged.append(df_curr)

        count += 1
        if count % 10000 == 0:
            logger.info("Iterate %i out of %i" % (count, n))

    df_merged.to_csv(output_merged_deduped_tsv_file, sep='\t', index=False)
    return df_merged
