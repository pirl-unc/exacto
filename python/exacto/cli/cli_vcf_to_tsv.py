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
and run Exacto 'vcf-to-tsv' command.
"""


import argparse
from ..main import *
from ..constants import *
from ..logging import get_logger


logger = get_logger(__name__)


def add_cli_vcf_to_tsv_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'vcf-to-tsv' parser.

    Parameters
    ----------
    sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns
    -------
    sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('vcf-to-tsv', help='Convert a VCF file to a TSV file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--vcf-file",
        dest="vcf_file",
        type=str,
        required=True,
        help="Input VCF file."
    )
    parser_required.add_argument(
        "--variant-calling-method",
        dest="variant_calling_method",
        type=str,
        required=True,
        choices=VariantCallingMethods.ALL,
        help="Variant calling method. "
             "Allowed options: %s."
             % (', '.join(VariantCallingMethods.ALL))
    )
    parser_required.add_argument(
        "--sequencing-platform",
        dest="sequencing_platform",
        type=str,
        required=True,
    )
    parser_required.add_argument(
        "--source-id",
        dest="source_id",
        type=str,
        required=True,
        help="Source ID (e.g. patient ID or cell line ID)."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser.set_defaults(which='vcf-to-tsv')
    return sub_parsers


def run_cli_vcf_to_tsv_from_parsed_args(args) -> None:
    """
    Run Exacto 'vcf-to-tsv' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                vcf_file
                variant_calling_method
                sequencing_platform
                output_tsv_file
                source_id
    """
    variants_list = run_exacto_vcf_to_tsv(
        vcf_file=args.vcf_file,
        variant_calling_method=args.variant_calling_method,
        sequencing_platform=args.sequencing_platform,
        source_id=args.source_id
    )
    df_variants = variants_list.to_dataframe()
    df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
