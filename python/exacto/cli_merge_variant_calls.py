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
and run Exacto 'merge' command.
"""


import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from .constants import *
from .default_parameters import *
from .logging import get_logger
from .main import *


logger = get_logger(__name__)


def add_cli_merge_variant_calls_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Adds 'merge-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('merge-variant-calls', help='Merge variant TSV files.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file",
        dest="tsv_file",
        action='append',
        required=True,
        help="Variants list TSV file."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser_required.add_argument(
        "--num-processes",
        dest="num_processes",
        type=int,
        default=MERGE_NUM_PROCESSES,
        required=True,
        help="Number of processes (default: %i)." % MERGE_NUM_PROCESSES
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--max-neighbor-distance",
        dest="max_neighbor_distance",
        type=int,
        required=False,
        default=MERGE_MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i)."
             % MERGE_MAX_NEIGHBOR_DISTANCE
    )
    parser.set_defaults(which='merge')
    return sub_parsers


def load_tsv_file_worker(tsv_file) -> VariantsList:
    """
    Loads a TSV file and returns a VariantsList object.

    Parameters
    ----------
    tsv_file        :   TSV file.

    Returns
    -------
    variants_list   :   VariantsList object.
    """
    variants_list = VariantsList.read_tsv_file(tsv_file=tsv_file)
    return variants_list


def run_cli_merge_variant_calls_from_parsed_args(args):
    """
    Run Exacto 'merge-variant-calls' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                'tsv_file'
                'output_tsv_file'
                'max_neighbor_distance'
                'num_processes'
    """
    # Step 1. Load variants lists
    logger.info("Started reading all TSV files")
    if args.num_processes > len(args.tsv_file):
        args.num_processes = len(args.tsv_file)
    tsv_files = np.array_split(args.tsv_file, args.num_processes)
    pool = mp.Pool(processes=args.num_processes)
    async_results = [pool.apply_async(load_tsv_file_worker, args=(tsv_file)) for tsv_file in tsv_files]
    pool.close()
    pool.join()
    variants_lists = [ar.get() for ar in async_results]
    logger.info("Finished reading all TSV files")

    # Step 2. Merge variants lists
    logger.info("Started merging all variant calls into one list")
    variants_list = run_exacto_merge_variant_calls(
        variants_lists=variants_lists,
        max_neighbor_distance=args.max_neighbor_distance
    )
    logger.info("Finished merging all variant calls into one list")

    # Step 3. Write to a TSV file
    variants_list.to_dataframe().to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
