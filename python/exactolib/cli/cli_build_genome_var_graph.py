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
and run Exacto 'build-genome-var-graph' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_build_genome_var_graph_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'build-genome-var-graph' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('build-genome-var-graph', help='Build a genome variation graph.')
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
        help="Reference genome FASTA file (variation graph backbone)."
    )
    parser_required.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=True,
        help="Output FASTA file."
    )
    parser_required.add_argument(
        "--sequence-prefix",
        dest="sequence_prefix",
        type=str,
        required=True,
        help="Sequence prefix."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--remove-unknown-bases",
        dest="remove_unknown_bases",
        type=str2bool,
        default=BUILD_GENOME_VAR_GRAPH_REMOVE_UNKNOWN_BASES,
        required=False,
        help="If 'yes', then unknown nucleotides ('N' or 'n') are removed (default: %s)."
             % BUILD_GENOME_VAR_GRAPH_REMOVE_UNKNOWN_BASES
    )
    parser_optional.add_argument(
        "--only-variant-sequences",
        dest="only_variant_sequences",
        type=str2bool,
        default=BUILD_GENOME_VAR_GRAPH_ONLY_VARIANT_SEQUENCES,
        required=False,
        help="If 'yes', then only variant sequences will be output (default: %s)."
             % BUILD_GENOME_VAR_GRAPH_ONLY_VARIANT_SEQUENCES
    )
    parser_optional.add_argument(
        "--graph-type",
        dest="graph_type",
        type=str,
        default=BUILD_GENOME_VAR_GRAPH_TYPE,
        required=False,
        help="Variation graph type. Either 'individual' or 'population' (default: %s)."
             % BUILD_GENOME_VAR_GRAPH_TYPE
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=BUILD_GENOME_VAR_GRAPH_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % BUILD_GENOME_VAR_GRAPH_NUM_THREADS
    )

    parser.set_defaults(which='build-genome-var-graph')
    return sub_parsers


def run_cli_build_genome_var_graph_from_parsed_args(args):
    """
    Run Exacto 'build-genome-var-graph' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    variants_tsv_file
                    fasta_file
                    output_fasta_file
                    sequence_prefix
                    num_threads
    """
    build_genome_variation_graph(
        df_variants=pl.read_csv(args.variants_tsv_file, separator='\t'),
        fasta_file=args.fasta_file,
        output_fasta_file=args.output_fasta_file,
        sequence_prefix=args.sequence_prefix,
        remove_unknown_bases=args.remove_unknown_bases,
        only_variant_sequences=args.only_variant_sequences,
        graph_type=args.graph_type,
        num_threads=args.num_threads,
        output_type=OutputType.FILE
    )
