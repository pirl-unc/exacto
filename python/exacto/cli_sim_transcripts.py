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
and run Exacto 'sim-transcripts' command.
"""


import argparse
import pandas as pd
from .constants import *
from .default_parameters import *
from .logging import get_logger
from .main import *


logger = get_logger(__name__)


def add_exacto_simulate_transcripts_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'sim-transcripts' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'sim-transcripts',
        help='Simulate transcripts.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        '--somatic_variants_tsv_file',
        dest='somatic_variants_tsv_file',
        type=str,
        required=True,
        help="Somatic (DNA and RNA) variants TSV file. "
             "Expected headers: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'variant_type', 'variant_sequence', 'clone_id' "
    )
    parser_required.add_argument(
        '--germline_variants_tsv_file',
        dest='germline_variants_tsv_file',
        type=str,
        required=True,
        help="Germline (DNA and RNA) variants TSV file. "
             "Expected headers: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'variant_type', 'variant_sequence', 'clone_id' "
    )
    parser_required.add_argument(
        '--gencode_reference_transcripts_gtf_file',
        dest='gencode_reference_transcripts_gtf_file',
        type=str,
        required=True,
        help="GENCODE reference transcripts GTF file."
    )
    parser_required.add_argument(
        '--output_transcripts_fasta_file',
        dest='output_transcripts_fasta_file',
        type=str,
        required=True,
        help="Output variants TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser.set_defaults(which='sim-transcripts')
    return sub_parsers


def run_exacto_simulate_transcripts_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'sim-transcripts' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                somatic_variants_tsv_file
                germline_variants_tsv_file
                gencode_reference_transcripts_gtf_file
                output_transcripts_fasta_file
    """
    # fasta = pysam.FastaFile(filename=args.fasta_file)
    # df_variants = run_exacto_simulate_variants(
    #     fasta=fasta
    # )
    # df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
    pass



