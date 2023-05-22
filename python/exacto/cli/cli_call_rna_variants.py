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
and run Exacto 'call-rna-variants' command.
"""


import argparse
import csv
import pandas as pd
import pysam
from ..constants import *
from ..default_parameters import *
from ..main import run_exacto_call_rna_variants


def add_cli_call_rna_variants_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'call-rna-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'call-rna-variants',
        help='Call RNA variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file."
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
        "--num-processes",
        dest="num_processes",
        type=int,
        default=CALL_VARIANTS_NUM_PROCESSES,
        required=False,
        help="Number of processes (default: %i)."
             % CALL_VARIANTS_NUM_PROCESSES
    )
    parser.set_defaults(which='call-rna-variants')
    return sub_parsers


def run_cli_call_rna_variants_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'call-rna-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                bam_file
                output_tsv_file
                num_processes
    """
    bam = pysam.AlignmentFile(args.bam_file)
    df_variants = run_exacto_call_rna_variants(
        bam=bam,
        num_processes=args.num_processes
    )
    df_variants.to_csv(
        args.output_tsv_file,
        sep='\t',
        quoting=csv.QUOTE_NONE,
        index=False
    )

