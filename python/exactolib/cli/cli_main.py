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
The purpose of this python3 script is to implement the primary Exacto command.
"""


import argparse
import exactolib
from .cli_identify import *
from ..logging import get_logger


logger = get_logger(__name__)


def init_arg_parser():
    """
    Initialize the input argument parser.

    Returns
    -------
    argparse.ArgumentParser object
    argparse.ArgumentParser subparsers object
    """
    arg_parser = argparse.ArgumentParser(
        description="Exacto: EXtracting Altered and Chimeric Transcripts Operations."
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(exactolib.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='Exacto sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_cli_identify_arg_parser(sub_parsers=sub_parsers)         # identify
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'identify':
        run_cli_identify_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)
