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
from .exacto_refine import *
from .exacto_annotate import *
from .exacto_merge import *
from .exacto_simulate_variants import *


logger = get_logger(__name__)


def init_arg_parser():
    """
    Initializes the input argument parser.

    Returns
    -------
    An instance of argparse.ArgumentParser
    An instance of argparse.ArgumentParser subparsers
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
    sub_parsers = add_exacto_simulate_variants_arg_parser(sub_parsers=sub_parsers)  # simulate_variants
    sub_parsers = add_exacto_refine_arg_parser(sub_parsers=sub_parsers)             # refine
    sub_parsers = add_exacto_annotate_arg_parser(sub_parsers=sub_parsers)           # annotate
    sub_parsers = add_exacto_merge_arg_parser(sub_parsers=sub_parsers)              # merge
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'call':
        a = 1
    elif args.which == 'graph':
        a = 1
    elif args.which == 'quantify':
        a = 1
    elif args.which == 'simulate_variants':
        run_exacto_simulate_variants_from_parsed_args(args=args)
    elif args.which == 'convert':
        a = 1
    elif args.which == 'refine':
        run_exacto_refine_from_parsed_args(args=args)
    elif args.which == 'annotate':
        run_exacto_annotate_from_parsed_args(args=args)
    elif args.which == 'merge':
        run_exacto_merge_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)

