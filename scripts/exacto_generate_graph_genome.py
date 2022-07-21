#!/usr/bin/python3

"""
The purpose of this python3 script is to generate a graph genome.

Last updated date: May 23, 2022

Author: Jin Seok (Andy) Lee
"""


import argparse
from exactolib.graph.genome import *


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="""
        Generates a graph genome based on a reference genome FASTA file
        and a list of somatic variants (TSV file).
        """
    )
    arg_parser.add_argument(
        "--reference_genome_fasta_file",
        dest="reference_genome_fasta_file",
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    arg_parser.add_argument(
        "--somatic_variants_tsv_file",
        dest="somatic_variants_tsv_file",
        type=str,
        required=True,
        help="Somatic variants TSV file."
    )
    arg_parser.add_argument(
        "--chromosomes",
        dest="chromosomes",
        type=str,
        required=True,
        help="Chromosomes (e.g. 'chr1 chr2 ... chrY')."
    )
    arg_parser.add_argument(
        "--output_graph_genome_file",
        dest="output_file",
        type=str,
        required=True,
        help="Output path (e.g. /<path>/output_directory/."
    )
    args = arg_parser.parse_args()
    args.chromosomes = args.chromosomes.split(' ')
    return (args)


if __name__ == '__main__':
    args = parse_args()


