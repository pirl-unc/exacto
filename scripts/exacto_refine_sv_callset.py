#!/usr/bin/python3

"""
The purpose of this python3 script is to refine variants.

Last updated date: July 15, 2022

Author: Jin Seok (Andy) Lee
"""


import argparse
from exactolib.logging import get_logger
from exactolib.constants import *
from exactolib.variant_refinement.structural_variants import *


logger = get_logger(__name__)


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="Refines a structural variant VCF file."
    )
    arg_parser._action_groups.pop()
    required = arg_parser.add_argument_group('required arguments')
    required.add_argument(
        "--vcf_file",
        dest="vcf_file",
        type=str,
        required=True,
        help="Input VCF file including path (e.g. /<path>/sample.vcf)."
    )
    required.add_argument(
        "--sv_calling_method",
        dest="sv_calling_method",
        type=str,
        required=True,
        choices=Constants.StructuralVariantCallingMethods.ALL,
        help="Structural variant calling method."
    )
    required.add_argument(
        "--sequencing_platform",
        dest="sequencing_platform",
        type=str,
        required=True,
        choices=Constants.SequencingPlatforms.ALL,
        help="Sequencing platform."
    )
    required.add_argument(
        "--blacklisted_regions_tsv_file",
        dest="blacklisted_regions_tsv_file",
        type=str,
        required=True,
        help="Blacklisted regions TSV file. Expected headers: 'chrom', 'chromStart', 'chromEnd'"
    )
    required.add_argument(
        "--gap_padding",
        dest="gap_padding",
        type=int,
        required=True,
        default=1E6,
        help="Number of bases to pad the gap regions (default: 1,000,000)."
    )
    required.add_argument(
        "--filter_values_to_include",
        dest="filter_values_to_include",
        action="append",
        nargs='+',
        required=True,
        help="Filter values to include (recommended: 'PASS')."
    )
    required.add_argument(
        "--min_total_coverage",
        dest="min_total_coverage",
        type=int,
        required=True,
        default=7,
        help="Minimum total coverage (default: 7)."
    )
    required.add_argument(
        "--min_variant_reads_count",
        dest="min_variant_reads_count",
        type=int,
        required=True,
        default=3,
        help="Minimum variant reads count (default: 3)."
    )
    required.add_argument(
        "--keep_only_precise",
        dest="keep_only_precise",
        type=bool,
        required=True,
        default=True,
        help="If true, only retains 'precise' variants (default: True)."
    )
    required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output refined TSV file (e.g. /<path>/output_directory/)."
    )
    required.add_argument(
        "--chromosomes",
        dest="chromosomes",
        type=str,
        required=True,
        action="append",
        nargs='+',
        help="Chromosomes to keep (e.g. --chromosomes chr1 chr2 chr3)."
    )
    return arg_parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.sv_calling_method == Constants.StructuralVariantCallingMethods.SNIFFLES2:
        df = refine_sniffles2_sv_callset(
            vcf_file=args.vcf_file,
            platform=args.sequencing_platform,
            blacklisted_regions_tsv_file=args.blacklisted_regions_tsv_file,
            filter_values_to_include=args.filter_values_to_include[0],
            min_total_coverage=args.min_total_coverage,
            min_variant_reads_count=args.min_variant_reads_count,
            keep_only_precise=args.keep_only_precise,
            gap_padding=args.gap_padding,
            chromosomes_to_keep=args.chromosomes[0]
        )
    if args.sv_calling_method == Constants.StructuralVariantCallingMethods.CUTESV:
        df = refine_cutesv_sv_callset(
            vcf_file=args.vcf_file,
            platform=args.sequencing_platform,
            blacklisted_regions_tsv_file=args.blacklisted_regions_tsv_file,
            filter_values_to_include=args.filter_values_to_include[0],
            min_total_coverage=args.min_total_coverage,
            min_variant_reads_count=args.min_variant_reads_count,
            keep_only_precise=args.keep_only_precise,
            gap_padding=args.gap_padding,
            chromosomes_to_keep=args.chromosomes[0]
        )
    if args.sv_calling_method == Constants.StructuralVariantCallingMethods.SVIM:
        df = refine_svim_sv_callset(
            vcf_file=args.vcf_file,
            platform=args.sequencing_platform,
            blacklisted_regions_tsv_file=args.blacklisted_regions_tsv_file,
            filter_values_to_include=args.filter_values_to_include[0],
            min_total_coverage=args.min_total_coverage,
            min_variant_reads_count=args.min_variant_reads_count,
            gap_padding=args.gap_padding,
            chromosomes_to_keep=args.chromosomes[0]
        )
    # if args.vcf_calling_method == Constants.VcfCallingMethods.PBSV:
    #     df = refine_pbsv_sv_callset(
    #         vcf_file=args.vcf_file,
    #         method=args.platform + '_' + args.vcf_calling_method,
    #         ucsc_gap_table_txt_file=args.ucsc_gap_table_txt_file,
    #         filter_values_to_include=args.filter_values_to_include[0],
    #         min_total_coverage=args.min_total_coverage,
    #         min_variant_reads_count=args.min_variant_reads_count,
    #         gap_padding=args.gap_padding,
    #         chromosomes_to_keep=args.chromosomes[0]
    #     )
    # if args.vcf_calling_method == Constants.VcfCallingMethods.DELLY2:
    #     df = refine_delly2_sv_callset(
    #         vcf_file=args.vcf_file,
    #         method=args.platform + '_' + args.vcf_calling_method,
    #         ucsc_gap_table_txt_file=args.ucsc_gap_table_txt_file,
    #         filter_values_to_include=args.filter_values_to_include[0],
    #         min_total_coverage=args.min_total_coverage,
    #         min_variant_reads_count=args.min_variant_reads_count,
    #         is_precise=args.is_precise,
    #         gap_padding=args.gap_padding,
    #         chromosomes_to_keep=args.chromosomes[0]
    #     )
    # if args.vcf_calling_method == Constants.VcfCallingMethods.LUMPY:
    #     df = refine_lumpy_callset(
    #         vcf_file=args.vcf_file,
    #         method=args.platform + '_' + args.vcf_calling_method,
    #         tumor_bam_file=args.tumor_bam_file,
    #         ucsc_gap_table_txt_file=args.ucsc_gap_table_txt_file,
    #         filter_values_to_include=args.filter_values_to_include[0],
    #         min_total_coverage=args.min_total_coverage,
    #         min_variant_reads_count=args.min_variant_reads_count,
    #         is_precise=args.is_precise,
    #         gap_padding=args.gap_padding,
    #         chromosomes_to_keep=args.chromosomes[0]
    #     )
    df.to_csv(args.output_tsv_file, sep='\t', index=False)
