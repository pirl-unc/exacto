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
and run Exacto 'sim-meiosis' command.
"""


import pandas as pd
from exacto.constants import *
from exacto.default_parameters import *


def add_exacto_simulate_meiosis_arg_parser(sub_parsers):
    """
    Adds 'simulate_meiosis' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'sim-meiosis',
        help='Simulates meiosis.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--germline_variants_tsv_file",
        dest="germline_variants_tsv_file",
        type=str,
        required=True,
        help="Paternal germline variants TSV file. "
             "The expected headers are: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'ref', 'alt', 'sv_type', 'genotype' "
             "(SNP, DEL, INS, INV, DUP, BND or TRA). "
             "If the 'sv_type' is 'SNP', leave 'chr_2' and 'pos_2' blank."
    )
    parser_required.add_argument(
        "--sex",
        dest="sex",
        type=str,
        choices=Sexes.ALL,
        required=True,
        help="Sex of simulated individual. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in Sexes.ALL))
    )
    parser_required.add_argument(
        "--num_meiototic_divisions",
        dest="num_meiototic_divisions",
        type=int,
        default=NUM_MEITOTIC_DIVISIONS,
        required=True,
        help="Number of meitotic divisions before the final gametes are sampled (default: %i). "
             "Each cell in each division generates 4 daughter cells (i.e. 4^n)."
             % NUM_MEITOTIC_DIVISIONS
    )
    parser_required.add_argument(
        "--num_sample_gametes",
        dest="num_sample_gametes",
        type=int,
        default=NUM_SAMPLE_GAMETES,
        required=True,
        help="Number of gametes to sample among the last generation of gametes (default: %i)."
             % NUM_SAMPLE_GAMETES
    )
    parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser.set_defaults(which='sim-meiosis')
    return sub_parsers


def run_exacto_simulate_meiosis_from_parsed_args(args):
    """
    Run Exacto 'simulate_meiosis' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                nucleic_acid_type
                fasta_file
                output_tsv_file
                num_snv
                num_insertion
                num_deletion
    """
    pass
    # fasta = pysam.FastaFile(filename=args.fasta_file)
    # df_variants = run_exacto_simulate_variants(
    #     fasta=fasta
    # )
    # df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

