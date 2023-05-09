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
The purpose of this python3 script is to implement the primary Exacto command.
"""


import argparse
import exacto
from .logging import get_logger
from .cli_annotate import *
from .cli_call_variants import *
from .cli_filter_variants import *
from .cli_merge_variant_calls import *
from .cli_sim_reads import *
from .cli_convert_vcf_to_tsv import *


logger = get_logger(__name__)


def init_arg_parser():
    """
    Initialize the input argument parser.

    Returns
    -------
    argparse.ArgumentParser object
    argparse.ArgumentParser subparsers object
    """
    arg_parser = argparse.ArgumentParser(
        description="Exacto: EXtracting And Counting Transcripts in Oncology."
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(exacto.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='Exacto sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_cli_convert_vcf_to_tsv_arg_parser(sub_parsers=sub_parsers)    # convert-vcf-to-tsv
    sub_parsers = add_cli_merge_variant_calls_arg_parser(sub_parsers=sub_parsers)   # merge-variant-calls
    sub_parsers = add_cli_filter_variants_arg_parser(sub_parsers=sub_parsers)       # filter-variants
    sub_parsers = add_cli_annotate_arg_parser(sub_parsers=sub_parsers)              # annotate
    sub_parsers = add_cli_call_variants_arg_parser(sub_parsers=sub_parsers)         # call-variants
    sub_parsers = add_cli_simulate_reads_arg_parser(sub_parsers=sub_parsers)        # sim-reads
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'convert-vcf-to-tsv':
        run_cli_convert_vcf_to_tsv_from_parsed_args(args=args)
    elif args.which == 'merge-variant-calls':
        run_cli_merge_variant_calls_from_parsed_args(args=args)
    elif args.which == 'filter-variants':
        run_cli_filter_variants_from_parsed_args(args=args)
    elif args.which == 'annotate':
        run_cli_annotate_from_parsed_args(args=args)
    elif args.which == 'call-variants':
        run_cli_call_variants_from_parsed_args(args=args)
    elif args.which == 'sim-reads':
        run_cli_sim_reads_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)
