import argparse
import pandas as pd
from exactolib.exactolibrs import format_gencode_gtf_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="Format GTF file for Exacto RNA variant calling."
    )

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "-i",
        dest="gtf_file",
        type=str,
        required=True,
        help="GTF file."
    )
    parser_required.add_argument(
        "-o",
        dest="output_txt_file",
        type=str,
        required=True,
        help="Output TXT file."
    )
    parser_required.add_argument(
        "-t",
        dest="num_threads",
        type=int,
        default=2,
        required=True,
        help="Number of threads."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "-s",
        dest="source",
        type=str,
        required=False,
        default='gencode',
        help="Reference FASTA source (default: 'gencode')."
    )
    return parser.parse_args()


def run():
    # Step 1. Read the input arguments
    args = parse_args()

    # Step 2. Format
    if args.source == 'gencode':
        pass
        # format_gencode_gtf_file(
        #     args.gtf_file,
        #     args.output_txt_file,
        #     args.num_threads
        # )
    else:
        raise Exception('Unsupported GTF file source: %s' % args.source)
