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
and run Exacto 'filter-variants' command.
"""


import argparse
from ..genomic_ranges_list import GenomicRangesList
from ..constants import VariantFilterSampleTypes
from ..default import *
from ..logging import get_logger
from ..main import run_exacto_filter_variants
from ..variant_filter import VariantFilter
from ..variants_list import VariantsList


logger = get_logger(__name__)


def add_cli_filter_variants_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Adds 'filter-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('filter-variants', help='Filter variants.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file",
        dest="tsv_file",
        type=str,
        required=True,
        help="Input variants TSV file."
    )
    parser_required.add_argument(
        "--case-sample-id",
        dest="case_sample_id",
        type=str,
        action='append',
        required=True,
        help="Case sample ID. Specify this parameter multiple times if there are "
             "multiple case sample IDs in the supplied TSV file. "
             "(e.g. --case-sample-id case_01 --case-sample-id case_02)."
    )
    parser_required.add_argument(
        "--output-filtered-tsv-file",
        dest="output_filtered_tsv_file",
        type=str,
        required=True,
        help="Output (filtered) TSV file."
    )
    parser_required.add_argument(
        "--output-rejected-tsv-file",
        dest="output_rejected_tsv_file",
        type=str,
        required=True,
        help="Output (rejected) TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        '--control-sample-id',
        dest='control_sample_id',
        type=str,
        action='append',
        required=False,
        help="Control sample ID. Specify this parameter multiple times if there are "
             "multiple control sample IDs in the supplied TSV file. "
             "(e.g. --control-sample-id control_01 --control-sample-id control_02)."
    )
    parser_optional.add_argument(
        '--filter',
        dest='filter',
        type=str,
        action='append',
        required=False,
        help='Variant filter conditions: '
             '"{case,control} {all,average,median,min,max,any} {attribute} {<,<=,>,>=,==,in} {value}". '
             'Example 1: "case all alternate_allele_read_count >= 3". '
             'Example 2: "case all chr_1 in ["chr1","chr2","chr3"]". '
             'Please refer to the Exacto documentation on how the filter semantics work.'
    )
    parser_optional.add_argument(
        "--excluded-regions-tsv-file",
        dest="excluded_regions_tsv_file",
        type=str,
        required=False,
        help="TSV files of regions to exclude. "
             "Variants with any variant calls where breakpoints are near the regions in this file will be removed. "
             "Expected headers: 'chromosome', 'start', 'end'."
    )
    parser_optional.add_argument(
        "--excluded-regions-padding",
        dest="excluded_regions_padding",
        type=int,
        required=False,
        default=FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING,
        help="Number of bases to pad each region in '--excluded-regions-tsv-files' (default: %i)."
             % FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING
    )
    parser_optional.add_argument(
        "--excluded-variants-tsv-file",
        dest="excluded_variants_tsv_file",
        type=str,
        required=False,
        help="TSV file of variants to exclude. Variants with any variant calls where "
             "breakpoints are near the variants in this file will be removed. Expected headers: "
             "'chromosome_1', 'position_1', 'chromosome_2', 'position_2', 'variant_type'. "
             "This parameter can be used to filter out germline variants."
    )
    parser_optional.add_argument(
        "--excluded-variants-padding",
        dest="excluded_variants_padding",
        type=int,
        required=False,
        default=FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING,
        help="Number of bases to pad the breakpoints of variants in '--excluded-variants-tsv-file' "
             "(default: %i)."
             % FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=FILTER_VARIANTS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)." % FILTER_VARIANTS_NUM_THREADS
    )
    parser.set_defaults(which='filter-variants')
    return sub_parsers


def run_cli_filter_variants_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'filter-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                tsv_file
                case_sample_id
                output_filtered_tsv_file
                output_rejected_tsv_file
                control_sample_id
                filter
                excluded_regions_tsv_file
                excluded_regions_padding
                excluded_variants_tsv_file
                excluded_variants_padding
                num_threads
    """
    # Step 1. Load variants list
    logger.info('Loading target VariantsList...')
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Load variant filters
    variant_filters = []
    if args.filter is not None:
        for curr_filter in args.filter:
            curr_filter = curr_filter.split(' ')
            if curr_filter[0] == VariantFilterSampleTypes.CASE:
                sample_ids = args.case_sample_id
            elif curr_filter[0] == VariantFilterSampleTypes.CONTROL:
                sample_ids = args.control_sample_id
            else:
                raise Exception('Unknown sample type for variant filter: %s.' % curr_filter[0])
            variant_filter = VariantFilter(
                quantifier=curr_filter[1],
                attribute=curr_filter[2],
                operator=curr_filter[3],
                value=curr_filter[4],
                sample_ids=sample_ids
            )
            variant_filters.append(variant_filter)
    logger.info('Loaded %i variant filters.' % len(variant_filters))

    # Step 3. Load excluded variants list
    if args.excluded_variants_tsv_file is not None:
        logger.info('Loading excluded VariantsList...')
        excluded_variants_list = VariantsList.read_tsv_file(tsv_file=args.excluded_variants_tsv_file)
    else:
        excluded_variants_list = None

    # Step 4. Load excluded regions list
    if args.excluded_regions_tsv_file is not None:
        excluded_regions_list = GenomicRangesList.read_tsv_file(tsv_file=args.excluded_regions_tsv_file)
        logger.info('Loaded %i excluded regions.' % excluded_regions_list.num_genomic_regions)
    else:
        excluded_regions_list = None

    # Step 4. Apply variant filters
    variants_list_filtered, variants_list_rejected = run_exacto_filter_variants(
        variants_list=variants_list,
        variant_filters=variant_filters,
        excluded_variants_list=excluded_variants_list,
        excluded_regions_list=excluded_regions_list,
        excluded_regions_padding=args.excluded_regions_padding,
        excluded_variants_padding=args.excluded_variants_padding,
        num_threads=args.num_threads
    )

    # Step 5. Write to TSV files
    logger.info('Started writing filtered and rejected variants')
    variants_list_filtered.to_dataframe().to_csv(args.output_filtered_tsv_file, sep='\t', index=False)
    variants_list_rejected.to_dataframe().to_csv(args.output_rejected_tsv_file, sep='\t', index=False)
    logger.info('Finished writing filtered and rejected variants')
