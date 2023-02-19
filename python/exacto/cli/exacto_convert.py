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
and run Exacto 'convert' command.
"""


import argparse
from ..constants import *
from ..default_parameters import *
from ..logging import get_logger
from ..main import *
from ..variants.vcf import *


logger = get_logger(__name__)


def add_exacto_convert_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'convert' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('convert', help='Convert a VCF file to a TSV file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        '--variant_class',
        type=str,
        required=True,
        choices=VariantClasses.ALL,
        help="Variant type (%s). "
             "If the input VCF file is of structural variants, specify '%s'. "
             "If the input VCF file is of SNVs and INDELs, specify '%s'."
             % (', '.join(f"'{item}'" for item in VariantClasses.ALL),
                VariantClasses.SV,
                VariantClasses.SNV_INDEL)
    )
    parser_required.add_argument(
        "--vcf_file",
        dest="vcf_file",
        type=str,
        required=True,
        help="Input VCF file."
    )
    parser_required.add_argument(
        "--variant_calling_method",
        dest="variant_calling_method",
        type=str,
        required=True,
        choices=VariantCallingMethods.ALL,
        help="Variant calling method. "
             "Supported options for 'sv': %s. "
             "Supported options for 'snv_indel': %s"
             % (', '.join(VariantCallingMethods.StructuralVariantCallingMethods.ALL),
                ', '.join(VariantCallingMethods.SmallVariantCallingMethods.ALL))
    )
    parser_required.add_argument(
        "--sequencing_platform",
        dest="sequencing_platform",
        type=str,
        required=True,
        choices=SequencingPlatforms.ALL,
        help="Sequencing platform."
    )
    parser_required.add_argument(
        "--sample_id",
        dest="sample_id",
        type=str,
        required=True,
        help="Sample ID."
    )
    parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output (refined) TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--tumor_sample_id",
        dest="tumor_sample_id",
        type=str,
        required=False,
        help="Tumor sample ID. "
             "This parameter must be specified if --variant_class is '%s'."
             % VariantClasses.SNV_INDEL
    )
    parser_optional.add_argument(
        "--normal_sample_id",
        dest="normal_sample_id",
        type=str,
        required=False,
        default='',
        help="Normal sample ID. This parameter must be specified if "
             "--variant_class is '%s' and variant calling was performed "
             "using a tumor and matched normal."
            % VariantClasses.SNV_INDEL
    )
    parser.set_defaults(which='convert')
    return sub_parsers


def run_exacto_convert_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'convert' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                variant_class
                vcf_file
                variant_calling_method
                sequencing_platform
                output_tsv_file
                sample_id
                tumor_sample_id
                normal_sample_id
    """
    # Step 1. Convert
    if args.variant_class == VariantClasses.SV:
        # Load VCF file
        if args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2:
            df_structural_variants = convert_sniffles2_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.CUTESV:
            df_structural_variants = convert_cutesv_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SVIM:
            df_structural_variants = convert_svim_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.PBSV:
            df_structural_variants = convert_pbsv_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id
            )
        else:
            raise Exception(
                "Invalid value for '--variant_calling_method': %s. "
                "Allowed '--variant_calling_method' values are %s "
                % (args.variant_calling_method,
                   ', '.join(f"'{item}'" for item in VariantCallingMethods.StructuralVariantCallingMethods.ALL))
            )

        # Write variants to a TSV file.
        df_structural_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
    elif args.variant_class == VariantClasses.SNV_INDEL:
        # Load VCF file
        if args.variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2:
            df_variants = convert_gatk4_mutect2_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id,
                tumor_sample_id=args.tumor_sample_id,
                normal_sample_id=args.normal_sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.STRELKA2_GERMLINE:
            df_variants = convert_strelka2_germline_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id,
                tumor_sample_id=args.tumor_sample_id,
                normal_sample_id=args.normal_sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT:
            df_variants = convert_deepvariant_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                sample_id=args.sample_id
            )
        else:
            raise Exception(
                "Invalid value for '--variant_calling_method': %s. "
                "Allowed '--variant_calling_method' values are %s "
                % (args.variant_calling_method,
                   ', '.join(f"'{item}'" for item in VariantCallingMethods.SmallVariantCallingMethods.ALL))
            )

        # Write variants to a TSV file.
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
    else:
        raise Exception(
            "Invalid value for '--variant_class': %s. "
            "Allowed '--variant_class' values are %s "
            % (args.variant_type,
               ', '.join(f"'{item}'" for item in VariantClasses.ALL))
        )
