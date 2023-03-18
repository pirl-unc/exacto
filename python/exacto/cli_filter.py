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
and run Exacto 'filter' command.
"""


import argparse
from .main import *


logger = get_logger(__name__)


def add_cli_filter_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'filter' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('filter', help='Filters variants.')
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
        '--filter',
        dest='filter',
        type=str,
        action='append',
        required=True,
        help='Filter conditions ("{dna,rna} {all,average,median,min,max,any} {attribute} {<,<=,>,>=,==,in} {value}").'
             'Example 1: "all alt_tumor_reads >= 3". '
             'Example 2: "all chr_1 in [chr1,chr2,chr3]". '
             'Please refer to the Exacto documentation on how the filter semantics work.'
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output (refined) TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--excluded-regions-tsv-files",
        dest="excluded_regions_tsv_files",
        type=str,
        required=False,
        help="TSV files of regions to exclude. "
             "Variant calls with breakpoints near the regions in this file will be removed. "
             "Expected headers: 'chrom', 'chromStart', 'chromEnd'."
    )
    parser_optional.add_argument(
        "--excluded-region-padding",
        dest="excluded_region_padding",
        type=int,
        required=False,
        default=EXCLUDED_REGION_PADDING,
        help="Number of bases to pad each region in '--excluded-regions-tsv-files' (default: %i)." % EXCLUDED_REGION_PADDING
    )
    parser_optional.add_argument(
        "--excluded-variants-tsv-files",
        dest="excluded_variants_tsv_files",
        type=str,
        required=False,
        nargs='+',
        help="TSV files of variants calls to explicitly exclude. "
             "Variant calls in '--tsv-file' that are close to a variant call in "
             "these files will be removed. Expected headers: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'variant_type'. "
             "This parameter can be used to filter out germline variants."
    )
    parser_optional.add_argument(
        "--excluded-variant-padding",
        dest="excluded_variant_padding",
        type=int,
        required=False,
        default=EXCLUDED_VARIANT_PADDING,
        help="Number of bases to pad the positions of variant calls in '--excluded-variants-tsv-files' "
             "(default: %i)."
             % EXCLUDED_VARIANT_PADDING
    )
    parser_optional.add_argument(
        "--enforce-variant-type-matching",
        dest="enforce_variant_type_matching",
        type=bool,
        required=False,
        default=ENFORCE_VARIANT_TYPE_MATCHING,
        help="If true, then the 'variant_type' of a variant call in "
             "'--tsv-file' must match a variant call in '--excluded-variants-tsv-files' "
             "for the variant call to be removed. If false, only the positions are considered "
             "when excluding a variant call (default: %r)."
    )
    parser.set_defaults(which='filter')
    return sub_parsers


def run_cli_filter_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'filter' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                tsv_file
                filter
                output_tsv_file
                excluded_regions_tsv_files
                excluded_region_padding
                excluded_variants_tsv_files
                excluded_variant_padding
                enforce_variant_type_matching
    """
    # Step 1. Load variant filters
    variant_filters = []
    for filter in args.filter:
        filter = filter.split(' ')
        variant_filter = VariantFilter(
            quantifier=filter[0],
            attribute=filter[1],
            operator=filter[2],
            value=filter[3]
        )
        variant_filters.append(variant_filter)

    # Step 2. Load excluded regions
    if args.excluded_regions_tsv_files is not None:
        df_excluded_regions = pd.DataFrame()
        for curr_tsv_file in args.excluded_regions_tsv_files:
            df_temp = pd.read_csv(curr_tsv_file, sep='\t')
            df_excluded_regions = pd.concat([df_excluded_regions, df_temp], axis=0)
    else:
        df_excluded_regions = pd.DataFrame()

    # Step 3. Load excluded variants
    if args.excluded_variants_tsv_files is not None:
        df_excluded_variants = pd.DataFrame()
        for tsv_file in args.excluded_variants_tsv_files:
            df_temp = pd.read_csv(tsv_file, sep='\t')
            df_excluded_variants = pd.concat([df_excluded_variants, df_temp], axis=0)
    else:
        df_excluded_variants = pd.DataFrame()

    # Step 4. Perform filtering
    variants_list = run_exacto_filter(
        variants_list=args.variants_list,
        df_excluded_variants=df_excluded_variants,
        df_excluded_regions=df_excluded_regions,
        variant_filters=variant_filters,
        excluded_region_padding=args.excluded_region_padding,
        excluded_variant_padding=args.excluded_variant_padding,
        enforce_variant_type_checking=args.enforce_variant_type_checking
    )

    # Step 5. Write to a TSV file
    variants_list.to_dataframe().to_csv(args.output_tsv_file, sep='\t', index=False)
