#!/usr/bin/python3

"""
The purpose of this python3 script is to merge SV TSV files.

Last updated date: July 20, 2022

Author: Jin Seok (Andy) Lee
"""


import argparse
from exactolib.utilities.merge_sv_tsv_files import *


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""Merges structural variant TSV files."""
    )
    arg_parser.add_argument(
        "--tsv_files",
        dest="tsv_files",
        action="append",
        nargs='+',
        required=True,
        help="List of SV TSV files including paths."
    )
    arg_parser.add_argument(
        "--output_merged_tsv_file",
        dest="output_merged_tsv_file",
        type=str,
        required=True,
        help="Output TSV file including path."
    )
    arg_parser.add_argument(
        "--output_merged_deduped_tsv_file",
        dest="output_merged_deduped_tsv_file",
        type=str,
        required=True,
        help="Output TSV file including path."
    )
    arg_parser.add_argument(
        "--methods_priority_list",
        dest="methods_priority_list",
        action="append",
        nargs='+',
        required=True,
        help="Priority list of methods."
    )
    arg_parser.add_argument(
        "--max_cluster_distance",
        dest="max_cluster_distance",
        type=int,
        required=True,
        default=10,
        help="""Maximum clustering distance (default: 10)."""
    )
    args = arg_parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    merge_sv_tsv_files(
        tsv_files=args.tsv_files[0],
        output_merged_tsv_file=args.output_merged_tsv_file,
        output_merged_deduped_tsv_file=args.output_merged_deduped_tsv_file,
        methods_priority_list=args.methods_priority_list[0],
        max_cluster_distance=args.max_cluster_distance
    )
