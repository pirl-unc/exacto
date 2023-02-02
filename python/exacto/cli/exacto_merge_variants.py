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


import pandas as pd
from ..logging import get_logger
from ..constants import *
from ..default_parameters import *
from ..main import *


logger = get_logger(__name__)


def add_exacto_merge_variants_arg_parser(sub_parsers):
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
        '--variant_class',
        type=str,
        required=True,
        choices=VariantClasses.ALL,
        help="Variant class (%s). "
             "If the input VCF file is of structural variants, specify '%s'. "
             "If the input VCF file is of SNVs and INDELs, specify '%s'."
             % (', '.join(f"'{item}'" for item in VariantClasses.ALL),
                VariantClasses.SV,
                VariantClasses.SNV_INDEL)
    )
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
    parser_required.add_argument(
        "--output_merged_deduped_tsv_file",
        dest="output_merged_deduped_tsv_file",
        type=str,
        required=True,
        help="Output merged and deduped TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--enforce_variant_type_matching",
        dest="enforce_variant_type_matching",
        type=bool,
        default=True,
        required=False,
        help="If true, variant type (i.e. 'sv_type' for structural variants and 'variant_type' for small variants) "
             "must match for 2 variants to be merged into one (default: True)."
    )
    parser_optional.add_argument(
        "--max_clustering_distance",
        dest="max_clustering_distance",
        type=int,
        required=False,
        default=MAX_SV_CLUSTER_DISTANCE,
        help="Maximum clustering distance (default: %i)."
             % MAX_SV_CLUSTER_DISTANCE
    )
    parser.set_defaults(which='merge-variants')
    return sub_parsers


def run_exacto_merge_variants_from_parsed_args(args):
    """
    Run Exacto 'merge-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                'variant_class'
                'tsv_files'
                'output_merged_tsv_file'
                'output_merged_deduped_tsv_file'
                'enforce_sv_type_matching'
                'max_clustering_distance'
    """
    if args.variant_class == VariantClasses.SV:
        list_df = []
        for curr_tsv_file in args.tsv_files:
            df_temp = pd.read_csv(curr_tsv_file, sep='\t')
            list_df.append(df_temp)
        df_merged, df_merged_deduped = run_exacto_merge_genomic_structural_variants(
            list_df=list_df,
            enforce_variant_type_matching=args.enforce_variant_type_matching,
            max_clustering_distance=args.max_clustering_distance
        )
        df_merged.to_csv(
            args.output_merged_tsv_file, sep='\t', index=False
        )
        df_merged_deduped.to_csv(args.output_merged_deduped_tsv_file,
                                 sep='\t',
                                 index=False)
    elif args.variant_class == VariantClasses.SNV_INDEL:
        list_df = []
        for curr_tsv_file in args.tsv_files:
            df_temp = pd.read_csv(curr_tsv_file, sep='\t')
            list_df.append(df_temp)
        df_merged, df_merged_deduped = run_exacto_merge_genomic_small_variants(
            list_df=list_df,
            enforce_variant_type_matching=args.enforce_variant_type_matching,
            max_clustering_distance=args.max_clustering_distance
        )
        df_merged.to_csv(
            args.output_merged_tsv_file, sep='\t', index=False
        )
        df_merged_deduped.to_csv(args.output_merged_deduped_tsv_file,
                                 sep='\t',
                                 index=False)
    else:
        raise Exception(
            "Invalid value for '--variant_class': %s. "
            "Allowed '--variant_class' values are %s "
            % (args.variant_class,
               ', '.join(f"'{item}'" for item in VariantClasses.ALL))
        )