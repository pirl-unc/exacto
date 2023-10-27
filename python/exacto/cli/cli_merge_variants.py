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
and run Exacto 'merge-variants' command.
"""


import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
import logging
from ..constants import *
from ..default import *
from ..logging import get_logger
from ..main import *


logger = get_logger(__name__)


def add_cli_merge_variants_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Adds 'merge-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('merge-variants', help='Merge variants.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file",
        dest="tsv_file",
        action='append',
        required=True,
        help="Variants list TSV file. "
             "This TSV file must follow Exacto's TSV format for this command to work properly."
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
        default=MERGE_VARIANTS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % MERGE_VARIANTS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--max-neighbor-distance",
        dest="max_neighbor_distance",
        type=int,
        required=False,
        default=MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE,
        help="Maximum neighbor distance (default: %i)."
             % MERGE_VARIANTS_MAX_NEIGHBOR_DISTANCE
    )
    parser.set_defaults(which='merge-variants')
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


def run_cli_merge_variants_from_parsed_args(args):
    """
    Run Exacto 'merge-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                tsv_file
                output_tsv_file
                num_threads
                max_neighbor_distance
    """
    # Step 1. Load variants lists
    logger.info("Started reading all TSV files")
    pool = mp.Pool(processes=args.num_threads)
    async_results = []
    for tsv_file in args.tsv_file:
        async_results.append(pool.apply_async(load_tsv_file_worker, args=(tsv_file,)))
    pool.close()
    pool.join()
    variants_lists = [async_result.get() for async_result in async_results]
    logger.info("Finished reading all TSV files")

    # Step 2. Merge variants lists
    logger.info("Started merging all variants into one list")
    variants_list = run_exacto_merge_variants(
        variants_lists=variants_lists,
        num_threads=args.num_threads,
        max_neighbor_distance=args.max_neighbor_distance
    )
    logger.info("Finished merging all variants into one list")

    # Step 3. Write to a TSV file
    df_variants_list = variants_list.to_dataframe()
    df_variants_list.sort_values(['variant_id'], inplace=True)
    df_variants_list.to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
