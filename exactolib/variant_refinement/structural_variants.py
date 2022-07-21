#!/usr/bin/python3

"""
The purpose of this python3 script is to implement functions that are used
to refine structural variants.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import os
import pandas as pd
import subprocess as sp
from exactolib.logging import get_logger
from exactolib.constants import *
from exactolib.utilities.vcf_to_dataframe import *


logger = get_logger(__name__)


def exclude_sv_from_bed_file(df, exclude_bed_file, breakpoint_padding=10):
    """
    Iterate through the dataframe df and
    excludes SV calls present in exclude_bed_file.

    Parameters
    ----------
    df                  : DataFrame with the following columns:
                          Chr_1, Pos_1, Chr_2, Pos_2, SV_Type
    exclude_bed_file    : Full path of BED file.
    breakpoint_padding  : Number of bases to pad for each breakpoint.

    Returns DataFrame
    """
    df_exclude = pd.read_csv(exclude_bed_file, sep='\t', header=None)
    df_exclude.columns = ['Chrom', 'Start', 'End', 'SV_Type']

    keep_ids = []
    exclude_ids = []
    for index, row in df.iterrows():
        curr_chr_1 = str(row['Chr_1'])
        curr_pos_1 = int(row['Pos_1'])
        curr_chr_2 = str(row['Chr_2'])
        curr_pos_2 = int(row['Pos_2'])
        curr_sv_type = str(row['SV_Type'])

        df_exclude_match = df_exclude.loc[
            (df_exclude['SV_Type'] == curr_sv_type) &
            (df_exclude['Chrom'] == curr_chr_1) & (df_exclude['Chrom'] == curr_chr_2) &
            (((df_exclude['Start'] >= (curr_pos_1 - breakpoint_padding)) &
              (df_exclude['Start'] <= (curr_pos_1 + breakpoint_padding)) &
              (df_exclude['End'] >= (curr_pos_2 - breakpoint_padding)) &
              (df_exclude['End'] <= (curr_pos_2 + breakpoint_padding))) |
              (df_exclude['Start'] >= (curr_pos_2 - breakpoint_padding)) &
              (df_exclude['Start'] <= (curr_pos_2 + breakpoint_padding)) &
              (df_exclude['End'] >= (curr_pos_1 - breakpoint_padding)) &
              (df_exclude['End'] <= (curr_pos_1 + breakpoint_padding)))
        ,:]

        if len(df_exclude_match) == 0:
            keep_ids.append(row['ID'])
        else:
            exclude_ids.append(row['ID'])
    df_refined = df.loc[df.ID.isin(keep_ids),:]
    df_excluded = df.loc[df.ID.isin(exclude_ids),:]
    return df_refined, df_excluded


def refine_sniffles2_sv_callset(vcf_file,
                                platform,
                                blacklisted_regions_tsv_file,
                                chromosomes_to_keep,
                                filter_values_to_include=['PASS'],
                                min_total_coverage=7,
                                min_variant_reads_count=3,
                                keep_only_precise=True,
                                gap_padding=1E6):
    """
    Refines a Sniffles VCF file and returns a dataframe of refined variants.

    Parameters
    ----------
    vcf_file                        : Full path of VCF file.
    platform                        : Sequencing platform.
    blacklisted_regions_tsv_file    : Full path of blacklisted regions TSV file.
                                      The expected column names are 'chrom', 'chromStart', 'chromEnd'.
    chromosomes_to_keep             : List of chromosomes to keep.
    filter_values_to_include        : List of FILTER values to include (default: ['PASS'])
    min_total_coverage              : Minimum total coverage (default: 7).
    min_variant_reads_count         : Minimum number of variants (support) reads (default: 3).
    keep_only_precise               : Retains PRECISE variants if true (default: True).
    gap_padding                     : Number of bases to pad for a variant to be
                                      considered in a gap region (default: 1E6).

    Returns
    -------
    DataFrame
    """
    # Step 1. Convert VCF file to DataFrame.
    df = convert_sniffles2_vcf_to_dataframe(
        vcf_file=vcf_file,
        method=platform + "_" + Constants.StructuralVariantCallingMethods.SNIFFLES2
    )
    logger.info('%i variants before refinement.' % len(df))

    # Step 2. Apply filters.
    df = df.loc[df['filter'].isin(filter_values_to_include),:] # filter
    df = df.loc[df['total_coverage'] >= min_total_coverage,:] # total coverage
    df = df.loc[df['variant_reads_count'] >= min_variant_reads_count,:] # variant reads count
    df = df.loc[df['chr_1'].isin(chromosomes_to_keep) & df['chr_2'].isin(chromosomes_to_keep),:]
    if keep_only_precise:
        df = df.loc[df['is_precise'] == True,:]

    # Step 3. Filter out variants where at least one of the two breakpoints
    # lies in the (padded) gap region.
    df_gaps = pd.read_csv(blacklisted_regions_tsv_file, sep='\t')
    df_gaps['start'] = df_gaps.apply(lambda row: int(row.chromStart - gap_padding), axis=1)
    df_gaps['end'] = df_gaps.apply(lambda row: int(row.chromEnd + gap_padding), axis=1)
    keep = []
    for index, row in df.iterrows():
        # Check if either breakpoint falls inside the (padded) gap region
        conditions = ((df_gaps['chrom'] == row['chr_1']) &
                      (df_gaps['start'] <= row['pos_1']) &
                      (df_gaps['end'] >= row['pos_1'])) | \
                     ((df_gaps['chrom'] == row['chr_2']) &
                      (df_gaps['start'] <= row['pos_2']) &
                      (df_gaps['end'] >= row['pos_2']))
        df_matched = df_gaps.loc[conditions,:]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df = df.loc[keep,:]
    logger.info('%i variants after refinement.' % len(df))
    return df


def refine_cutesv_sv_callset(vcf_file,
                             platform,
                             blacklisted_regions_tsv_file,
                             chromosomes_to_keep,
                             filter_values_to_include=['PASS'],
                             min_total_coverage=7,
                             min_variant_reads_count=3,
                             keep_only_precise=True,
                             gap_padding=1E6):
    """
    Refines a cuteSV VCF file and returns a dataframe of refined variants.

    Parameters
    ----------
    vcf_file                        : Full path of VCF file.
    platform                        : Sequencing platform.
    blacklisted_regions_tsv_file    : Full path of blacklisted regions TSV file.
                                      The expected column names are 'chrom', 'chromStart', 'chromEnd'.
    chromosomes_to_keep             : List of chromosomes to keep.
    filter_values_to_include        : List of FILTER values to include (default: ['PASS'])
    min_total_coverage              : Minimum total coverage (default: 7).
    min_variant_reads_count         : Minimum number of variants (support) reads (default: 3).
    keep_only_precise               : Retains PRECISE variants if true (default: True).
    gap_padding                     : Number of bases to pad for a variant to be
                                      considered in a gap region (default: 1E6).

    Returns
    -------
    DataFrame
    """
    # Step 1. Convert VCF file to DataFrame.
    df = convert_cutesv_vcf_to_dataframe(
        vcf_file=vcf_file,
        method=platform + "_" + Constants.StructuralVariantCallingMethods.CUTESV
    )
    logger.info('%i variants before refinement.' % len(df))

    # Step 2. Apply filters.
    df = df.loc[df['filter'].isin(filter_values_to_include),:] # filter
    df = df.loc[df['total_coverage'] >= min_total_coverage,:] # total coverage
    df = df.loc[df['variant_reads_count'] >= min_variant_reads_count,:] # variant reads count
    df = df.loc[df['chr_1'].isin(chromosomes_to_keep) & df['chr_2'].isin(chromosomes_to_keep),:]
    if keep_only_precise:
        df = df.loc[df['is_precise'] == True,:]

    # Step 3. Filter out variants where at least one of the two breakpoints
    # lies in the (padded) gap region.
    df_gaps = pd.read_csv(blacklisted_regions_tsv_file, sep='\t')
    df_gaps['start'] = df_gaps.apply(lambda row: int(row.chromStart - gap_padding), axis=1)
    df_gaps['end'] = df_gaps.apply(lambda row: int(row.chromEnd + gap_padding), axis=1)
    keep = []
    for index, row in df.iterrows():
        # Check if either breakpoint falls inside the (padded) gap region
        conditions = ((df_gaps['chrom'] == row['chr_1']) &
                      (df_gaps['start'] <= row['pos_1']) &
                      (df_gaps['end'] >= row['pos_1'])) | \
                     ((df_gaps['chrom'] == row['chr_2']) &
                      (df_gaps['start'] <= row['pos_2']) &
                      (df_gaps['end'] >= row['pos_2']))
        df_matched = df_gaps.loc[conditions,:]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df = df.loc[keep,:]
    logger.info('%i variants after refinement.' % len(df))
    return df


def refine_svim_sv_callset(vcf_file,
                           platform,
                           blacklisted_regions_tsv_file,
                           chromosomes_to_keep,
                           filter_values_to_include=['PASS'],
                           min_total_coverage=7,
                           min_variant_reads_count=3,
                           gap_padding=1E6):
    """
    Refines a SVIM VCF file and returns a dataframe of refined variants.

    Parameters
    ----------
    vcf_file                        : Full path of VCF file.
    platform                        : Sequencing platform.
    blacklisted_regions_tsv_file    : Full path of blacklisted regions TSV file.
                                      The expected column names are 'chrom', 'chromStart', 'chromEnd'.
    chromosomes_to_keep             : List of chromosomes to keep.
    filter_values_to_include        : List of FILTER values to include (default: ['PASS'])
    min_total_coverage              : Minimum total coverage (default: 7).
    min_variant_reads_count         : Minimum number of variants (support) reads (default: 3).
    gap_padding                     : Number of bases to pad for a variant to be
                                      considered in a gap region (default: 1E6).

    Returns
    -------
    DataFrame
    """
    # Step 1. Convert VCF file to DataFrame.
    df = convert_svim_vcf_to_dataframe(
        vcf_file=vcf_file,
        method=platform + "_" + Constants.StructuralVariantCallingMethods.SVIM
    )
    logger.info('%i variants before refinement.' % len(df))

    # Step 2. Apply filters.
    df = df.loc[df['filter'].isin(filter_values_to_include),:] # filter
    df = df.loc[df['total_coverage'] >= min_total_coverage,:] # total coverage
    df = df.loc[df['variant_reads_count'] >= min_variant_reads_count,:] # variant reads count
    df = df.loc[df['chr_1'].isin(chromosomes_to_keep) & df['chr_2'].isin(chromosomes_to_keep),:]

    # Step 3. Filter out variants where at least one of the two breakpoints
    # lies in the (padded) gap region.
    df_gaps = pd.read_csv(blacklisted_regions_tsv_file, sep='\t')
    df_gaps['start'] = df_gaps.apply(lambda row: int(row.chromStart - gap_padding), axis=1)
    df_gaps['end'] = df_gaps.apply(lambda row: int(row.chromEnd + gap_padding), axis=1)
    keep = []
    for index, row in df.iterrows():
        # Check if either breakpoint falls inside the (padded) gap region
        conditions = ((df_gaps['chrom'] == row['chr_1']) &
                      (df_gaps['start'] <= row['pos_1']) &
                      (df_gaps['end'] >= row['pos_1'])) | \
                     ((df_gaps['chrom'] == row['chr_2']) &
                      (df_gaps['start'] <= row['pos_2']) &
                      (df_gaps['end'] >= row['pos_2']))
        df_matched = df_gaps.loc[conditions,:]
        if len(df_matched) == 0:
            keep.append(True)
        else:
            keep.append(False)
    df = df.loc[keep,:]
    logger.info('%i variants after refinement.' % len(df))
    return df
