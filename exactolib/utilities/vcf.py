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
handling VCF files.
"""


import vcf
import pandas as pd
from collections import defaultdict
from ..logging import get_logger
from ..constants import *


logger = get_logger(__name__)


def convert_gatk4_mutect2_vcf_to_dataframe(
        vcf_file: str,
        sequencing_platform: str,
        tumor_sample_id: str) -> pd.DataFrame:
    """
    Convert a GATK4-Mutect2 VCF file
    (called with tumor and normal samples) to a DataFrame.

    Parameters
    ----------
    vcf_file    :   VCF file.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chrom'
    'pos'
    'ref'
    'alt'
    'variant_type'
    'filter'
    'as_sb_table'
    'ecnt'
    'germq'
    'mbq'
    'mfrl'
    'mmq'
    'mpos'
    'nalod'
    'nlod'
    'popaf'
    'tlod'
    'tumor_reference_reads_count'
    'tumor_variant_reads_count'
    'tumor_total_coverage'
    'variant_allele_fraction'
    'normal_reference_reads_count'
    'normal_total_coverage'
    'tumor_genotype'
    'normal_genotype'
    """
    data = defaultdict(list)
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        chrom = str(record.CHROM)                                                           # chromosome
        pos = int(record.POS)                                                               # position
        if 'chr' not in chrom:                                                              # make sure 'chr' is in chrom
            chrom = 'chr' + chrom
        ref = str(record.REF)                                                               # reference allele
        alt = record.ALT                                                                    # alternate allele
        if len(alt) > 1:
            continue
        else:
            alt = str(record.ALT[0])
        if len(ref) == 1 and len(alt) == 1:                                                 # variant type
            variant_type = 'snv'
        elif len(ref) == 1 and len(alt) > 1:
            variant_type = 'insertion'
        elif len(ref) > 1 and len(alt) == 1:
            variant_type = 'deletion'
        elif len(ref) >= 2 and len(alt) >= 2:
            variant_type = 'mnv'
        else:
            logger.info("Unknown variant type. REF: %s. ALT: %s" % (ref, alt))
        if len(record.FILTER) == 0:                                                         # filter value
            filter = 'PASS'
        else:
            filter = record.FILTER[0]
        as_sb_table = record.INFO['AS_SB_TABLE']                                            # allele-specific forward/reverse read counts for strand bias tests.
        ecnt = record.INFO['ECNT']                                                          # number of events in this haplotype
        germq = record.INFO['GERMQ']                                                        # phred-scale quality that alt alleles are not germline variants
        mbq = record.INFO['MBQ'][0]                                                         # median base quality by allele
        mfrl = record.INFO['MFRL'][0]                                                       # median fragment length by allele
        mmq = record.INFO['MMQ'][0]                                                         # median mapping quality by allele
        mpos = record.INFO['MPOS'][0]                                                       # median distance from end of read
        nalod = record.INFO['NALOD'][0]                                                     # negative log 10 odds of artifact in normal with same allele fraction as tumor
        nlod = record.INFO['NLOD'][0]                                                       # normal log 10 odds of artifact in normal with same allele fraction as tumor
        popaf = record.INFO['POPAF'][0]                                                     # negative log 10 population allele frequencies of alt alleles
        tlod = record.INFO['TLOD'][0]                                                       # log 10 likelihood ratio score of variant existing versus not existing
        if record.samples[0].sample == tumor_sample_id:                                     # tumor and normal sample record index
            tumor_sample_record_idx = 0
            normal_sample_record_idx = 1
        else:
            tumor_sample_record_idx = 1
            normal_sample_record_idx = 0
        tumor_reference_reads_count = record.samples[tumor_sample_record_idx]['AD'][0]      # tumor reference reads count
        tumor_variant_reads_count = record.samples[tumor_sample_record_idx]['AD'][1]        # tumor variant reads count
        tumor_total_coverage = record.samples[tumor_sample_record_idx]['DP']                # tumor total coverage
        normal_reference_reads_count = record.samples[normal_sample_record_idx]['AD'][0]    # normal reference reads count
        normal_variant_reads_count = record.samples[normal_sample_record_idx]['AD'][1]      # normal variant reads count
        normal_total_coverage = record.samples[normal_sample_record_idx]['DP']              # normal total coverage
        if tumor_variant_reads_count == 0 or \
                tumor_total_coverage == 0 or \
                normal_variant_reads_count > 0:
            continue
        variant_allele_fraction = tumor_variant_reads_count / tumor_total_coverage          # variant allele fraction
        tumor_genotype = str(record.samples[tumor_sample_record_idx]['GT'])                 # tumor genotype
        normal_genotype = str(record.samples[normal_sample_record_idx]['GT'])               # normal genotype

        # Append to data
        data['id'].append(record.ID)
        data['variant_calling_method'].append(VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2)
        data['sequencing_platform'].append(sequencing_platform)
        data['chrom'].append(chrom)
        data['pos'].append(pos)
        data['ref'].append(ref)
        data['alt'].append(alt)
        data['variant_type'].append(variant_type)
        data['filter'].append(filter)
        data['as_sb_table'].append(as_sb_table)
        data['ecnt'].append(ecnt)
        data['germq'].append(germq)
        data['mbq'].append(mbq)
        data['mfrl'].append(mfrl)
        data['mmq'].append(mmq)
        data['mpos'].append(mpos)
        data['nalod'].append(nalod)
        data['nlod'].append(nlod)
        data['popaf'].append(popaf)
        data['tlod'].append(tlod)
        data['tumor_reference_reads_count'].append(tumor_reference_reads_count)
        data['tumor_variant_reads_count'].append(tumor_variant_reads_count)
        data['tumor_total_coverage'].append(tumor_total_coverage)
        data['variant_allele_fraction'].append(variant_allele_fraction)
        data['normal_reference_reads_count'].append(normal_reference_reads_count)
        data['normal_total_coverage'].append(normal_total_coverage)
        data['tumor_genotype'].append(tumor_genotype)
        data['normal_genotype'].append(normal_genotype)

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning dataframe' % len(df))
    return df


def convert_sniffles2_vcf_to_dataframe(vcf_file: str, sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a Sniffles2 VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
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
        try:                                                                # standard deviation of structural variant size
            sv_size_stdev = float(record.INFO['STDEV_LEN'])
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
        data['variant_calling_method'].append(VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2)
        data['sequencing_platform'].append(sequencing_platform)
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
        data['sv_size_stdev'].append(sv_size_stdev)
        data['variant_reads_count'].append(variant_reads_count)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        data['variant_allele_fraction'].append(vaf)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['sv_pos_stdev'].append(sv_pos_stdev)
        data['coverage'].append(coverage)
        data['query_alignment_length_adjusted_mismatches_mean_count'].append(nm)
        data['support_long'].append(support_long)
        data['ci_pos'].append("unknown")
        data['ci_len'].append("unknown")
        data['std_span'].append("unknown")
        data['tandem_duplication_copy_number'].append("unknown")
        data['strand_reads'].append("unknown")
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_cutesv_vcf_to_dataframe(vcf_file: str, sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a CuteSV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
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
        data['variant_calling_method'].append(VariantCallingMethods.StructuralVariantCallingMethods.CUTESV)
        data['sequencing_platform'].append(sequencing_platform)
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
        data['sv_size_stdev'].append("unknown")
        data['variant_reads_count'].append(variant_reads_count)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        data['variant_allele_fraction'].append(vaf)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['sv_pos_stdev'].append("unknown")
        data['coverage'].append("unknown")
        data['query_alignment_length_adjusted_mismatches_mean_count'].append("unknown")
        data['support_long'].append("unknown")
        data['ci_pos'].append(ci_pos)
        data['ci_len'].append(ci_len)
        data['std_span'].append("unknown")
        data['tandem_duplication_copy_number'].append("unknown")
        data['strand_reads'].append("unknown")
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_svim_vcf_to_dataframe(vcf_file: str, sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a SVIM VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'sv_type'
    'sv_size'
    'sv_pos_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
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
            sv_pos_stdev = str(record.INFO['STD_POS'])
            if sv_pos_stdev is None:
                sv_pos_stdev = "unknown"
        except:
            sv_pos_stdev = "unknown"
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
        data['variant_calling_method'].append(VariantCallingMethods.StructuralVariantCallingMethods.SVIM)
        data['sequencing_platform'].append(sequencing_platform)
        data['chr_1'].append(chr_1)
        data['pos_1'].append(pos_1)
        data['chr_2'].append(chr_2)
        data['pos_2'].append(pos_2)
        data['ref'].append(ref_allele)
        data['alt'].append(alt_allele)
        data['quality_score'].append(quality_score)
        data['filter'].append(filter)
        data['is_precise'].append("unknown")
        data['sv_type'].append(sv_type)
        data['sv_size'].append(sv_size)
        data['sv_size_stdev'].append("unknown")
        data['variant_reads_count'].append(variant_reads_count)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        data['variant_allele_fraction'].append(vaf)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['sv_pos_stdev'].append(sv_pos_stdev)
        data['coverage'].append("unknown")
        data['query_alignment_length_adjusted_mismatches_mean_count'].append("unknown")
        data['support_long'].append("unknown")
        data['ci_pos'].append("unknown")
        data['ci_len'].append("unknown")
        data['std_span'].append(std_span)
        data['tandem_duplication_copy_number'].append(tandem_duplication_copy_number)
        data['strand_reads'].append("unknown")
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df


def convert_pbsv_vcf_to_dataframe(vcf_file: str, sequencing_platform: str) -> pd.DataFrame:
    """
    Convert a PBSV VCF file to a DataFrame.

    Parameters
    ----------
    vcf_file                :   Path to VCF file.
    sequencing_platform     :   Sequencing platform.

    Returns
    -------
    DataFrame with the following columns:
    'id'
    'variant_calling_method'
    'sequencing_platform'
    'chr_1'
    'pos_1'
    'chr_2'
    'pos_2'
    'ref'
    'alt'
    'quality_score'
    'filter'
    'is_precise'
    'sv_type'
    'sv_size'
    'sv_size_stdev'
    'variant_reads_count'
    'reference_reads_count'
    'total_coverage'
    'variant_allele_fraction'
    'read_ids'
    'strand'
    'insertion_sequence'
    'genotype'
    'genotype_quality'
    'sv_pos_stdev'
    'coverage'
    'query_alignment_length_adjusted_mismatches_mean_count'
    'support_long'
    'ci_pos'
    'ci_len'
    'std_span'
    'tandem_duplication_copy_number'
    'strand_reads'
    """
    data = defaultdict(list)
    count = 0
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    included_mate_ids = set()

    for record in vcf_reader:
        chr_1 = str(record.CHROM)                                           # chromosome 1
        pos_1 = int(record.POS)                                             # position 1
        if str(record.INFO['SVTYPE']) == 'BND':                             # chromosome 2
            curr_id = str(record.ID).split("-")[1]
            curr_id = curr_id.split(":")[0]
            chr_2 = str(curr_id)
        else:
            chr_2 = chr_1
        if 'chr' not in chr_1:                                              # make sure 'chr' is in chr_1 and chr_2
            chr_1 = 'chr' + chr_1
        if 'chr' not in chr_2:
            chr_2 = 'chr' + chr_2
        if str(record.INFO['SVTYPE']) == 'BND':                             # position 2
            curr_id = str(record.ID).split("-")[1]
            curr_id = curr_id.split(":")[1]
            pos_2 = int(curr_id)
        else:
            pos_2 = int(record.INFO['END'])
        ref_allele = ""                                                     # reference allele
        alt_allele = ""                                                     # alternate allele
        quality_score = record.QUAL                                         # quality score
        if len(record.FILTER) == 0:                                         # filter
            filter = "PASS"
        else:
            filter = record.FILTER[0]
        if 'IMPRECISE' in record.INFO:                                      # is breakpoint precise?
            is_precise = False
        else:
            is_precise = True
        sv_type = str(record.INFO['SVTYPE'])                                # SV type
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
        variant_reads_count = int(record.samples[0]['AD'][0])               # number of supporting reads
        variant_read_ids = ''                                               # read IDs
        if sv_type == 'BND':                                                # strand
            alt_val = str(record.ALT[0])
            if (alt_val[0:2] == 'A[') or (alt_val[0:2] == 'C[') or (alt_val[0:2] == 'T[') or (alt_val[0:2] == 'G['):
                strand = '+-'
            elif (alt_val[0:2] == 'A]') or (alt_val[0:2] == 'C]') or (alt_val[0:2] == 'T]') or (alt_val[0:2] == 'G]'):
                strand = '-+'
            elif (alt_val[-2:] == ']A') or (alt_val[-2:] == ']C') or (alt_val[-2:] == ']T') or (alt_val[-2:] == ']G'):
                strand = '+-'
            elif (alt_val[-2:] == '[A') or (alt_val[-2:] == '[C') or (alt_val[-2:] == '[T') or (alt_val[-2:] == '[G'):
                strand = '-+'
            else:
                strand = 'unknown'
        else:
            strand = 'unknown'
        if sv_type == 'INS':                                                # insertion sequence
            insertion_sequence = str(record.ALT[0]).upper()
        else:
            insertion_sequence = ''
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
            strand_reads = str(record.samples[0]['SAC'])                    # forward and reverse strand reads in each allele
        except:
            strand_reads = ''
        if sv_type == 'BND':
            ci_pos = str(record.INFO['CIPOS'])
        else:
            ci_pos = ''

        if sv_type == 'BND':
            # Check if current ID has been included
            if str(record.INFO['MATEID'][0]) in included_mate_ids:
                continue
            included_mate_ids.add(record.ID)

        # Append to data
        data['id'].append(str(record.ID))
        data['variant_calling_method'].append(VariantCallingMethods.StructuralVariantCallingMethods.PBSV)
        data['sequencing_platform'].append(sequencing_platform)
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
        data['sv_size_stdev'].append("unknown")
        data['variant_reads_count'].append(variant_reads_count)
        data['reference_reads_count'].append(reference_reads_count)
        data['total_coverage'].append(total_coverage)
        data['variant_allele_fraction'].append(vaf)
        data['read_ids'].append(variant_read_ids)
        data['strand'].append(strand)
        data['insertion_sequence'].append(insertion_sequence)
        data['genotype'].append(genotype)
        data['genotype_quality'].append(genotype_quality)
        data['sv_pos_stdev'].append("unknown")
        data['ci_pos'].append("unknown")
        data['ci_len'].append("unknown")
        data['std_span'].append("unknown")
        data['tandem_duplication_copy_number'].append("unknown")
        data['strand_reads'].append(strand_reads)
        count += 1

    df = pd.DataFrame(data)
    logger.info('%i rows in the returning DataFrame.' % len(df))
    return df
