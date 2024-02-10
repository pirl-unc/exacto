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
and run Paradgm 'identify' command.
"""


import argparse
import pandas as pd
from collections import defaultdict
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import *


logger = get_logger(__name__)


def add_cli_identify_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'identify' parser.

    Parameters
    ----------
    sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns
    -------
    sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('identify', help='Identify RNA variants in a BAM file.')
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
        "--num-threads",
        dest="num_threads",
        type=int,
        default=IDENTIFY_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % IDENTIFY_NUM_THREADS
    )
    parser_optional.add_argument(
        '--chromosomes',
        dest='chromosomes',
        type=str,
        nargs='+',
        required=False,
        help='Chromosomes. If unspecified, Exacto identifies variants in all chromosomes.'
    )
    parser_optional.add_argument(
        "--min-reads",
        dest="min_reads",
        type=int,
        default=IDENTIFY_MIN_READS,
        required=False,
        help="Minimum number of supporting reads (default: %i)."
             % IDENTIFY_MIN_READS
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=IDENTIFY_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % IDENTIFY_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-ins-size-proportion",
        dest="min_ins_size_proportion",
        type=float,
        default=IDENTIFY_MIN_INS_SIZE_PROPORTION,
        required=False,
        help="Minimum insertion size proportion between two insertions (default: %f). "
             "Size proportion = smaller insertion size / longer insertion size."
             % IDENTIFY_MIN_INS_SIZE_PROPORTION
    )
    parser_optional.add_argument(
        "--max-ins-norm-edit-distance",
        dest="max_ins_norm_edit_distance",
        type=float,
        default=IDENTIFY_MAX_INS_NORM_EDIT_DISTANCE,
        required=False,
        help="Maximum insertion normalized edit (Levenshtein) distance (default: %f). "
             "Normalized edit distance = edit distance / longer insertion size."
             % IDENTIFY_MAX_INS_NORM_EDIT_DISTANCE
    )
    parser_optional.add_argument(
        "--min-del-size-proportion",
        dest="min_del_size_proportion",
        type=float,
        default=IDENTIFY_MIN_DEL_SIZE_PROPORTION,
        required=False,
        help="Minimum deletion size proportion between two deletions (default: %f). "
             "Size proportion = smaller deletion size / longer deletion size."
             % IDENTIFY_MIN_DEL_SIZE_PROPORTION
    )

    parser.set_defaults(which='identify')
    return sub_parsers


def run_cli_identify_from_parsed_args(args) -> None:
    """
    Run Exacto 'identify' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                bam_file
                output_tsv_file
                num_threads
                chromosome
                min_reads
                min_mapping_quality
                min_ins_size_proportion
                max_ins_norm_edit_distance
                min_del_size_proportion
    """
    variant_calls = identify_rna_variants(
        bam_file=args.bam_file,
        min_reads=args.min_reads,
        min_mapping_quality=args.min_mapping_quality,
        num_threads=args.num_threads,
        min_ins_size_proportion=args.min_ins_size_proportion,
        max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
        min_del_size_proportion=args.min_del_size_proportion,
        chromosomes=args.chromosomes
    )
    data = defaultdict(list)
    for variant_call in variant_calls:
        for key, value in variant_call.to_dict().items():
            data[key].append(value)
    df_variant_calls = pd.DataFrame(data)
    df_variant_calls.to_csv(args.output_tsv_file, sep='\t', index=False)
