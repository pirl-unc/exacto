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
and run Exacto 'diff-kmers' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_diff_kmers_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'diff-kmers' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('diff-kmers', help='Diff k-mers.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--query-fasta-file",
        dest="query_fasta_file",
        type=str,
        required=True,
        help="Query FASTA file."
    )
    parser_required.add_argument(
        "--reference-fasta-file",
        dest="reference_fasta_file",
        type=str,
        required=True,
        help="Reference peptides (e.g. GENCODE) FASTA file."
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
        "--min-k",
        dest="min_k",
        type=int,
        default=DIFF_KMERS_MIN_K,
        required=False,
        help="Minimum k-mer k size (default: %i)."
             % DIFF_KMERS_MIN_K
    )
    parser_optional.add_argument(
        "--max-k",
        dest="max_k",
        type=int,
        default=DIFF_KMERS_MAX_K,
        required=False,
        help="Maximum k-mer k size (default: %i)."
             % DIFF_KMERS_MAX_K
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_PEPTIDE_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_PEPTIDE_VARS_NUM_THREADS
    )

    parser.set_defaults(which='diff-kmers')
    return sub_parsers


def run_cli_diff_kmers_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'diff-kmers' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    query_fasta_file
                    reference_fasta_file
                    reference_fasta_file
                    output_tsv_file
                    output_fasta_file
                    min_k
                    max_k
                    num_threads
    """
    pass
    # # Read the query FASTA file peptide sequences
    # query_sequences = {} # key = peptide ID, value = peptide sequence
    # query_fasta = pysam.FastaFile(args.query_fasta_file)
    # for sequence_name in query_fasta.references:
    #     sequence = query_fasta.fetch(sequence_name)
    #     query_sequences[sequence_name] = sequence
    #
    # # Read the reference FASTA file peptide sequences
    # reference_sequences = {} # key = peptide ID, value = peptide sequence
    # reference_fasta = pysam.FastaFile(args.reference_fasta_file)
    # for sequence_name in reference_fasta.references:
    #     sequence = reference_fasta.fetch(sequence_name)
    #     reference_sequences[sequence_name] = sequence
    #
    # # Identify unique kmers
    # df_unique_kmers = diff_kmers(
    #     query_sequences=query_sequences,
    #     reference_sequences=reference_sequences,
    #     min_k=args.min_k,
    #     max_k=args.max_k
    # )
    #
    # # Write to TSV file
    # df_unique_kmers.to_csv(args.output_tsv_file, sep='\t', index=False)
    #
