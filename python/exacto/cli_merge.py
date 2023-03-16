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
The purpose of this python3 script is to create parser
and run Exacto 'merge' command.
"""


import argparse
import pandas as pd
from .constants import *
from .default_parameters import *
from .logging import get_logger
from .main import *


logger = get_logger(__name__)


def add_cli_merge_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Adds 'merge' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('merge', help='Merge variant TSV files.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-files",
        dest="tsv_files",
        nargs='+',
        required=True,
        help="List of variant TSV files."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--enforce_variant_type_matching",
        dest="enforce_variant_type_matching",
        type=bool,
        default=True,
        required=False,
        help="If true, variant type (i.e. 'sv_type' for structural variants and 'variant_type' for small variants) "
             "must match for 2 variants to be merged into one (default: True)."
    )
    parser_optional.add_argument(
        "--max-neighbor-distance",
        dest="max_neighbor_distance",
        type=int,
        required=False,
        default=MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i)."
             % MAX_NEIGHBOR_DISTANCE
    )
    parser.set_defaults(which='merge')
    return sub_parsers


def run_cli_merge_from_parsed_args(args):
    """
    Run Exacto 'merge' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                'tsv_files'
                'output_tsv_file'
                'enforce_variant_type_matching'
                'max_neighbor_distance'
    """
    variants_lists = []
    for tsv_file in args.tsv_files:
        variants_lists.append(VariantsList.read_tsv_file(tsv_file=tsv_file))
    variants_list = run_exacto_merge(
        variants_lists=variants_lists,
        enforce_variant_type_matching=args.enforce_variant_type_matching,
        max_neighbor_distance=args.max_variant_merge_distance
    )
    variants_list.to_dataframe().to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
