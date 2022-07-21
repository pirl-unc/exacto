#!/usr/bin/python3

"""
The purpose of this python3 script is to implement VCF class and functions
related to handling a VCF file.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""

import vcf
import pandas as pd
from collections import defaultdict
from exactolib.logging import get_logger
from exactolib.constants import *


logger = get_logger(__name__)


def convert_sniffles2_vcf_to_dataframe(vcf_file, method):
    """
    Convert a Sniffles2 VCF file to a DataFrame.

    Args
    ----
    vcf_file    :   Path to VCF file.
    method      :   Method.

    Returns
    -------
    DataFrame
    """
    data = defaultdict(list)
    count = 0
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        chr_1 = str(record.CHROM)                                           # chromosome 1
        pos_1 = int(record.POS)                                             # position 1
        if str(record.INFO['SVTYPE']) == 'BND':                             # chromosome 2
            chr_2 = str(record.INFO['CHR2'])
        else:
            chr_2 = chr_1
        if 'chr' not in chr_1:                                              # make sure 'chr' is in chr_1 and chr_2
            chr_1 = 'chr' + chr_1
        if 'chr' not in chr_2:
            chr_2 = 'chr' + chr_2
        if str(record.INFO['SVTYPE']) == 'BND':                             # position 2
            alt_val = str(record.ALT[0]).split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            pos_2 = int(alt_val)
        else:
            pos_2 = int(record.INFO['END'])
        ref_allele = str(record.REF).upper()                                # reference allele
        alt_allele = str(record.ALT[0]).upper()                             # alternate allele
        quality_score = float(record.QUAL)                                  # quality score
        if len(record.FILTER) == 0:                                         # filter
            filter = "PASS"
        else:
            filter = record.FILTER[0]
        if 'PRECISE' in record.INFO:                                        # is breakpoint precise?
            is_precise = True
        else:
            is_precise = False
        sv_type = str(record.INFO['SVTYPE'])                                # SV type
        if sv_type == 'BND':                                                # SV size
            if chr_1 == chr_2:
                sv_size = abs(pos_2 - pos_1)
            else:
                sv_size = -1
        else:
            sv_size = abs(int(record.INFO['SVLEN']))
        variant_reads_count = int(record.INFO['SUPPORT'])                   # number of supporting reads
        variant_read_ids = ','.join(record.INFO['RNAMES'])                  # read IDs
        coverage = ','.join(record.INFO['COVERAGE'])                        # coverage (upstream, start, center, end, downstream)
        strand = record.INFO['STRAND']                                      # strand
        if sv_type == 'INS':                                                # insertion sequence
            insertion_sequence = str(record.ALT[0]).upper()
        else:
            insertion_sequence = ''
        nm = int(record.INFO['NM'][0])                                      # mean number of query alignment length adjusted mismatches of supporting reads
        try:                                                                # standard deviation of structural variant length
            sv_len_stdev = float(record.INFO['STDEV_LEN'])
        except:
            sv_len_stdev = -1.0
        sv_pos_stdev = float(record.INFO['STDEV_POS'])                      # standard deviation of structural variant start position
        if sv_type == 'INS':                                                # number of soft-clipped reads putatively supporting the long insertion SV
            support_long = int(record.INFO['SUPPORT_LONG'])
        else:
            support_long = -1
        try:
            vaf = float(record.INFO['AF'])                                  # variant allele fraction
        except:
            vaf = -1.0
        genotype = record.samples[0]['GT']                                  # genotype
        genotype_quality = record.samples[0]['GQ']                          # genotype quality
        reference_reads_count = record.samples[0]['DR']                     # number of reference reads
        total_coverage = reference_reads_count + variant_reads_count  # total coverage

        # Append to data
        data['id'].append(str(record.ID))
        data['chr_1'].append(chr_1)
        data['pos_1'].append(pos_1)
        data['chr_2'].append(chr_2)
        data['pos_2'].append(pos_2)
        data['ref'].append(ref_allele)
        data['alt'].append(alt_allele)
        data['quality_score'].append(quality_score)
        data['filter'].append(filter)
        data['is_precise'].append(is_precise)
        data['sv_type'].append(sv_type)
        data['sv_size'].append(sv_size)
        data['method'].append(method)
        data['variant_reads_count'].append(variant_reads_count)
        data['read_ids'].append(variant_read_ids)
        data['coverage'].append(coverage)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['query_alignment_length_adjusted_mismatches_mean_count'].append(nm)
        data['sv_len_stdev'].append(sv_len_stdev)
        data['sv_pos_stdev'].append(sv_pos_stdev)
        data['support_long'].append(support_long)
        data['variant_allele_fraction'].append(vaf)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_cutesv_vcf_to_dataframe(vcf_file, method):
    """
    Convert a CuteSV VCF file to a DataFrame.

    Args
    ----
    vcf_file    :   Path to VCF file.
    method      :   Method.

    Returns
    -------
    DataFrame
    """
    data = defaultdict(list)
    count = 0
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        chr_1 = str(record.CHROM)                                           # chromosome 1
        pos_1 = int(record.POS)                                             # position 1
        if str(record.INFO['SVTYPE']) == 'BND':                             # chromosome 2
            alt_val = str(record.ALT[0]).split(":")[0]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            chr_2 = str(alt_val)
        else:
            chr_2 = chr_1
        if 'chr' not in chr_1:                                              # make sure 'chr' is in chr_1 and chr_2
            chr_1 = 'chr' + chr_1
        if 'chr' not in chr_2:
            chr_2 = 'chr' + chr_2
        if str(record.INFO['SVTYPE']) == 'BND':                             # position 2
            alt_val = str(record.ALT[0]).split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            pos_2 = int(alt_val)
        else:
            pos_2 = int(record.INFO['END'])
        ref_allele = str(record.REF).upper()                                # reference allele
        alt_allele = str(record.ALT[0]).upper()                             # alternate allele
        quality_score = record.QUAL                                         # quality score
        if len(record.FILTER) == 0:                                         # filter
            filter = "PASS"
        else:
            filter = record.FILTER[0]
        if 'PRECISE' in record.INFO:                                        # is breakpoint precise?
            is_precise = True
        else:
            is_precise = False
        sv_type = str(record.INFO['SVTYPE'])                                # SV type
        if sv_type == 'BND':                                                # SV size
            if chr_1 == chr_2:
                sv_size = abs(pos_2 - pos_1)
            else:
                sv_size = -1
        else:
            sv_size = abs(int(record.INFO['SVLEN']))
        variant_reads_count = int(record.INFO['RE'])                        # number of supporting reads
        variant_read_ids = ','.join(record.INFO['RNAMES'])                  # read IDs
        if sv_type == 'BND':                                                # strand
            alt_val = str(record.ALT[0])
            if alt_val[0:2] == 'N[':                                        # first part of strand
                strand = '+-'
            elif alt_val[0:2] == 'N]':
                strand = '-+'
            elif alt_val[-2:] == ']N':
                strand = '+-'
            elif alt_val[-2:] == '[N':
                strand = '-+'
            else:
                strand = 'unknown'
        else:
            try:
                strand = record.INFO['STRAND'][0]
            except:
                strand = 'unknown'
        if sv_type == 'INS':                                                # insertion sequence
            insertion_sequence = str(record.ALT[0])[1:].upper()
        else:
            insertion_sequence = ''
        try:                                                                # confidence interval around POS for impreicse variants
            ci_pos = ','.join(record.INFO['CIPOS'])
        except:
            ci_pos = "unknown"
        try:                                                                # confidence interval around inserted / deleted material between breakends
            ci_len = ','.join(record.INFO['CILEN'])
        except:
            ci_len = "unknown"
        genotype = record.samples[0]['GT']                                  # genotype
        genotype_quality = record.samples[0]['GQ']                          # genotype quality
        reference_reads_count = record.samples[0]['DR']                     # number of reference reads
        if reference_reads_count is None:
            reference_reads_count = -1
            total_coverage = -1
            vaf = -1.0
        else:
            total_coverage = reference_reads_count + variant_reads_count  # total coverage
            vaf = float(variant_reads_count) / float(total_coverage)

        # Append to data
        data['id'].append(str(record.ID))
        data['chr_1'].append(chr_1)
        data['pos_1'].append(pos_1)
        data['chr_2'].append(chr_2)
        data['pos_2'].append(pos_2)
        data['ref'].append(ref_allele)
        data['alt'].append(alt_allele)
        data['quality_score'].append(quality_score)
        data['filter'].append(filter)
        data['is_precise'].append(is_precise)
        data['sv_type'].append(sv_type)
        data['sv_size'].append(sv_size)
        data['method'].append(method)
        data['variant_reads_count'].append(variant_reads_count)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['ci_pos'].append(ci_pos)
        data['ci_len'].append(ci_len)
        data['variant_allele_fraction'].append(vaf)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_svim_vcf_to_dataframe(vcf_file, method):
    """
    Convert a SVIM VCF file to a DataFrame.

    Args
    ----
    vcf_file    :   Path to VCF file.
    method      :   Method.

    Returns
    -------
    DataFrame
    """
    data = defaultdict(list)
    count = 0
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        chr_1 = str(record.CHROM)                                           # chromosome 1
        pos_1 = int(record.POS)                                             # position 1
        if str(record.INFO['SVTYPE']) == 'BND':                             # chromosome 2
            alt_val = str(record.ALT[0]).split(":")[0]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            chr_2 = str(alt_val)
        else:
            chr_2 = chr_1
        if 'chr' not in chr_1:                                              # make sure 'chr' is in chr_1 and chr_2
            chr_1 = 'chr' + chr_1
        if 'chr' not in chr_2:
            chr_2 = 'chr' + chr_2
        if str(record.INFO['SVTYPE']) == 'BND':                             # position 2
            alt_val = str(record.ALT[0]).split(":")[1]
            alt_val = alt_val.replace("[", "")
            alt_val = alt_val.replace("]", "")
            alt_val = alt_val.replace("N", "")
            pos_2 = int(alt_val)
        else:
            pos_2 = int(record.INFO['END'])
        ref_allele = ""                                                     # reference allele
        alt_allele = ""                                                     # alternate allele
        quality_score = record.QUAL                                         # quality score
        if len(record.FILTER) == 0:                                         # filter
            filter = "PASS"
        else:
            filter = record.FILTER[0]
        is_precise = 'unknown'                                              # is breakpoint precise?
        sv_type = str(record.INFO['SVTYPE'])                                # SV type
        if "DUP" in sv_type:
            sv_type = 'DUP'
        if sv_type == 'BND':                                                # SV size
            if chr_1 == chr_2:
                sv_size = abs(pos_2 - pos_1)
            else:
                sv_size = -1
        else:
            try:
                sv_size = abs(int(record.INFO['SVLEN']))
            except:
                sv_size = -1
        variant_reads_count = int(record.INFO['SUPPORT'])                   # number of supporting reads
        variant_read_ids = ','.join(record.INFO['READS'])                   # read IDs
        if sv_type == 'BND':                                                # strand
            alt_val = str(record.ALT[0])
            if alt_val[0:2] == 'N[':                                        # first part of strand
                strand = '+-'
            elif alt_val[0:2] == 'N]':
                strand = '-+'
            elif alt_val[-2:] == ']N':
                strand = '+-'
            elif alt_val[-2:] == '[N':
                strand = '-+'
            else:
                strand = 'unknown'
        else:
            strand = 'unknown'
        if sv_type == 'INS':                                                # insertion sequence
            insertion_sequence = ','.join(record.INFO['SEQS'])
        else:
            insertion_sequence = ''
        try:                                                                # standard deviation in span of merged SV signatures
            std_span = str(record.INFO['STD_SPAN'])
            if std_span is None:
                std_span = "unknown"
        except:
            std_span = "unknown"
        try:                                                                # standard deviation in position of merged SV signatures
            std_pos = str(record.INFO['STD_POS'])
            if std_pos is None:
                std_pos = "unknown"
        except:
            std_pos = "unknown"
        genotype = record.samples[0]['GT']                                  # genotype
        genotype_quality = 'unknown'                                        # genotype quality
        total_coverage = record.samples[0]['DP']                            # total coverage
        if total_coverage is None:
            reference_reads_count = -1
            total_coverage = -1
            vaf = -1.0
        else:
            reference_reads_count = total_coverage - variant_reads_count
            vaf = float(variant_reads_count) / float(total_coverage)
        try:
            tandem_duplication_copy_number = int(record.samples[0]['CN'])   # copy number of tandem duplication (2 for one additional copy)
        except:
            tandem_duplication_copy_number = -1.0

        # Append to data
        data['id'].append(str(record.ID))
        data['chr_1'].append(chr_1)
        data['pos_1'].append(pos_1)
        data['chr_2'].append(chr_2)
        data['pos_2'].append(pos_2)
        data['ref'].append(ref_allele)
        data['alt'].append(alt_allele)
        data['quality_score'].append(quality_score)
        data['filter'].append(filter)
        data['is_precise'].append(is_precise)
        data['sv_type'].append(sv_type)
        data['sv_size'].append(sv_size)
        data['method'].append(method)
        data['variant_reads_count'].append(variant_reads_count)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['std_span'].append(std_span)
        data['std_pos'].append(std_pos)
        data['variant_allele_fraction'].append(vaf)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        data['tandem_duplication_copy_number'].append(tandem_duplication_copy_number)
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df

