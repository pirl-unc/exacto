#!/usr/bin/python3

"""
The purpose of this python3 script is to annotate
structural variants (TSV file) using pyensembl.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import argparse
from exactolib.logging import get_logger
from exactolib.variant_annotation.annotate_structural_variants import *


logger = get_logger(__name__)


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""Annotates structural variants (TSV file) using pyensembl."""
    )
    arg_parser.add_argument(
        "--ensembl_release",
        dest="ensembl_release",
        type=int,
        required=True,
        help="Ensembl release version number (recommended: 75 for GRCh37 106 for GRCh38)."
    )
    arg_parser.add_argument(
        "--tsv_file",
        dest="tsv_file",
        type=str,
        required=True,
        help="""
        Input TSV file including path (e.g. /<path>/sample.vcf). 
        The expected column headers are: 'chr_1', 'pos_1', 'chr_2', 'pos_2', 'sv_type' (DEL, INS, INV, DUP, BND).
        """
    )
    arg_parser.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file including full path."
    )

    args = arg_parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    annotate_using_pyensembl(
        tsv_file=args.tsv_file,
        ensembl_release=args.ensembl_release,
        output_tsv_file=args.output_tsv_file
    )
