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
and run Exacto 'translate-structs' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_translate_structs_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'translate-structs' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('translate-structs', help='Translate transcript structures to primary structures.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--transcript-structures-tsv-file",
        dest="transcript_structures_tsv_file",
        type=str,
        required=True,
        help="Input transcript structures TSV file."
    )
    parser_required.add_argument(
        "--rna-variant-calls-tsv-file",
        dest="rna_variant_calls_tsv_file",
        type=str,
        required=True,
        help="Input RNA variant calls TSV file."
    )
    parser_required.add_argument(
        "--integrated-variants-tsv-file",
        dest="integrated_variants_tsv_file",
        type=str,
        required=True,
        help="Input integrated variants TSV file."
    )
    parser_required.add_argument(
        "--strategy",
        dest="strategy",
        type=str,
        choices=[str(TranslationStrategy.LONGEST_ORF), str(TranslationStrategy.ALL_ORFS)],
        required=True,
        help="Translation strategy."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser_required.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=True,
        help="Output FASTA file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=TRANSLATE_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % TRANSLATE_NUM_THREADS
    )
    parser.set_defaults(which='translate-structs')
    return sub_parsers


def run_cli_translate_structs_from_parsed_args(args) -> None:
    """
    Run Exacto 'translate-structs' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    fastq_file
                    strategy
                    output_tsv_file
                    output_fasta_file
                    num_threads
                    gzip
    """
    translate_structures(
        transcript_structures_tsv_file=args.transcript_structures_tsv_file,
        rna_variant_calls_tsv_file=args.rna_variant_calls_tsv_file,
        integrated_variants_tsv_file=args.integrated_variants_tsv_file,
        strategy=TranslationStrategy(args.strategy),
        output_tsv_file=args.output_tsv_file,
        output_fasta_file=args.output_fasta_file,
        num_threads=args.num_threads
    )
