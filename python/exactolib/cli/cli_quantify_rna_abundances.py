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
and run Exacto 'quantify-rna-abundances' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_quantify_rna_abundances_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'quantify-rna-abundances' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('quantify-rna-abundances', help='Quantify RNA abundances.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--clusters-tsv-file",
        dest="clusters_tsv_file",
        type=str,
        required=True,
        help="Clusters TSV file. Expected columns: 'cluster_id', 'read_name'. "
             "This file is one of the outputs from exacto cluster-rna-reads."
    )
    parser_required.add_argument(
        "--cluster-reference-matches-tsv-file",
        dest="cluster_reference_matches_tsv_file",
        type=str,
        required=True,
        help="Cluster reference matches TSV file. Expected columns: 'cluster_id', 'reference_gene_id', 'reference_transcript_id'. "
             "This file is one of the outputs from exacto cluster-rna-reads."
    )
    parser_required.add_argument(
        "--cluster-splice-junctions-tsv-file",
        dest="cluster_splice_junctions_tsv_file",
        type=str,
        required=True,
        help="Cluster splice junctions TSV file. Expected columns: 'cluster_id', 'chromosome_1', 'chromosome_2', 'position_1', 'position_2', 'strand_1', 'strand_2'. "
             "This file is one of the outputs from exacto cluster-rna-reads."
    )
    parser_required.add_argument(
        "--cluster-variants-tsv-file",
        dest="cluster_variants_tsv_file",
        type=str,
        required=True,
        help="Cluster variants TSV file. Expected columns: 'cluster_id', 'chromosome_1', 'position_1', 'strand_1', 'operation_type_1', 'chromosome_2', 'position_2', 'strand_2', 'operation_type_2', 'sequence', 'variant_type'. "
             "This file is one of the outputs from exacto cluster-rna-reads."
    )
    parser_required.add_argument(
        "--consensus-rna-reference-matches-tsv-file",
        dest="consensus_rna_reference_matches_tsv_file",
        type=str,
        required=True,
        help="Consensus RNA reference transcript matches TSV file. Expected columns: 'assembled_transcript_name', 'reference_gene_name', 'reference_transcript_id'. "
             "This file is one of the outputs from exacto call-rna-transcript-vars."
    )
    parser_required.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        required=True,
        help="Output directory."
    )
    parser_required.add_argument(
        "--output-prefix",
        dest="output_prefix",
        type=str,
        required=True,
        help="Output prefix."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--pseudo-count",
        dest="pseudo_count",
        type=float,
        default=QUANTIFY_RNA_ABUNDANCES_PSEUDO_COUNT,
        required=False,
        help="Pseudo count added to each cluster's unique read count when the abundances are "
             "initialized (default: %s)." % QUANTIFY_RNA_ABUNDANCES_PSEUDO_COUNT
    )
    parser_optional.add_argument(
        "--max-iter",
        dest="max_iter",
        type=int,
        default=QUANTIFY_RNA_ABUNDANCES_MAX_ITER,
        required=False,
        help="Maximum number of expectation-maximization iterations (default: %i)."
             % QUANTIFY_RNA_ABUNDANCES_MAX_ITER
    )
    parser_optional.add_argument(
        "--tol",
        dest="tol",
        type=float,
        default=QUANTIFY_RNA_ABUNDANCES_TOL,
        required=False,
        help="Convergence tolerance: iteration stops once the largest per-cluster change in "
             "relative abundance falls below this (default: %s)." % QUANTIFY_RNA_ABUNDANCES_TOL
    )

    parser.set_defaults(which='quantify-rna-abundances')
    return sub_parsers


def run_cli_quantify_rna_abundances_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'quantify-rna-abundances' command using parameters from parsed arguments.
    """
    quantify_rna_abundances(
        clusters_tsv_file=args.clusters_tsv_file,
        cluster_reference_matches_tsv_file=args.cluster_reference_matches_tsv_file,
        cluster_splice_junctions_tsv_file=args.cluster_splice_junctions_tsv_file,
        cluster_variants_tsv_file=args.cluster_variants_tsv_file,
        consensus_rna_reference_matches_tsv_file=args.consensus_rna_reference_matches_tsv_file,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        pseudo_count=args.pseudo_count,
        max_iter=args.max_iter,
        tol=args.tol
    )
