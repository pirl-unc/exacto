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
and run Exacto 'translate-seqs' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_translate_seqs_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'translate-seqs' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('translate-seqs', help='Translate transcript sequences in a long-read RNA-seq FASTX file or RNA sequence to peptide sequences.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')

    # Mutually exclusive group for input
    input_group = parser_required.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--fastx-file",
        dest="fastx_file",
        type=str,
        help="Input FASTA or FASTQ file."
    )
    input_group.add_argument(
        "--sequence",
        dest="sequence",
        type=str,
        help="Directly input a transcript sequence."
    )

    parser_required.add_argument(
        "--strategy",
        dest="strategy",
        type=str,
        choices=[str(TranslationStrategy.LONGEST_ORF), str(TranslationStrategy.ALL_ORFS)],
        required=True,
        help="Translation strategy."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--start-codons",
        dest="start_codons",
        type=str,
        nargs='+',
        default=["AUG"],
        required=False,
        help="One or more start codons (default: AUG). Pass multiple as: --start-codons AUG GUG CUG."
    )
    parser_optional.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=False,
        help="Output TSV file."
    )
    parser_optional.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=False,
        help="Output FASTA file."
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=TRANSLATE_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % TRANSLATE_NUM_THREADS
    )

    parser.set_defaults(which='translate-seqs')
    return sub_parsers


def run_cli_translate_seqs_from_parsed_args(args) -> None:
    """
    Run Exacto 'translate-seqs' command using parameters from parsed arguments.
    """
    if args.fastx_file:
        if args.output_tsv_file is None:
            raise Exception('--output-tsv-file must be specified.')
        if args.output_fasta_file is None:
            raise Exception('--output-fasta-file must be specified.')
        logger.info("%i reads in total in the FASTX file." % count_reads_in_fastx_file(fastx_file=args.fastx_file))

        translate_fastx_file(
            fastx_file=args.fastx_file,
            output_fasta_file=args.output_fasta_file,
            output_tsv_file=args.output_tsv_file,
            strategy=TranslationStrategy(args.strategy),
            start_codons=args.start_codons,
            num_threads=args.num_threads
        )
    elif args.sequence:
        translations = translate_sequence(
            rna_sequence=args.sequence,
            strategy=args.strategy,
            start_codons=args.start_codons,
        )
        logger.info('Translated peptide sequence(s):')
        logger.info('[orf_start:orf_end] [sequence]')
        for (sequence, orf_start, orf_end) in translations:
            logger.info('%i:%i %s' % (orf_start, orf_end, sequence))
    else:
        raise Exception("Unexpected error: one of the following should have been specified: --fastx-file or --sequence.")
