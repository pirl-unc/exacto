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
and run Exacto 'build-var-graph' command.
"""


import argparse
from ..constants import *
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_build_var_graph_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'build-var-graph' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('build-var-graph', help='Build a variation graph.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--variants-tsv-file",
        dest="variants_tsv_file",
        type=str,
        required=True,
        help="Variants TSV file. Expected columns: 'variant_id', 'chromosome_1', 'position_1', 'operation_1', 'strand_1', 'chromosome_2', 'position_2', 'operation_2', 'strand_2', 'sequence'.",
    )
    parser_required.add_argument(
        "--fasta-file",
        dest="fasta_file",
        type=str,
        required=True,
        help="Fasta file (backbone of variation graph)."
    )
    parser_required.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=True,
        help="Output fasta file."
    )

    # Optional arguments
    parser.set_defaults(which='build-var-graph')
    return sub_parsers


def run_cli_build_var_graph_from_parsed_args(args):
    """
    Run Exacto 'build-var-graph' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    variants_tsv_file
                    fasta_file
                    output_fasta_file
    """
    pass
    # build_variation_graph(
    #     variants_tsv_file=args.variants_tsv_file,
    #     fasta_file=args.fasta_file,
    #     output_fasta_file=args.output_fasta_file
    # )
