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
import pandas as pd
from ..constants import *
from ..default_parameters import *
from ..logging import get_logger
from ..main import *


logger = get_logger(__name__)


def add_exacto_merge_annotations_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'merge-annotations' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('merge-annotations', help='Merge annotations.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv_files",
        dest="tsv_files",
        nargs='+',
        required=True,
        help="List of TSV files. "
             "Expected columns in each TSV file: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'sv_type'."
    )
    parser_required.add_argument(
        "--output_merged_tsv_file",
        dest="output_merged_tsv_file",
        type=str,
        required=True,
        help="Output merged TSV file."
    )

    parser.set_defaults(which='merge-annotations')
    return sub_parsers


def run_exacto_merge_annotations_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'merge-annotations' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                'tsv_files'
                'output_merged_tsv_file'
    """
    list_df = []
    for curr_tsv_file in args.tsv_files:
        df_temp = pd.read_csv(curr_tsv_file, sep='\t')
        list_df.append(df_temp)
    df_merged = run_exacto_merge_annotations(list_df=list_df)
    df_merged.to_csv(args.output_merged_tsv_file, sep='\t', index=False)
