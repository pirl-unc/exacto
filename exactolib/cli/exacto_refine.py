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
and run Exacto 'refine' command.
"""


from ..default_parameters import *
from ..constants import *
from ..main import *
from ..logging import get_logger
from ..utilities.vcf_utils import *
from ..variant_refinement.structural_variants import *
from ..variant_refinement.small_variants import *


logger = get_logger(__name__)


def add_exacto_refine_arg_parser(sub_parsers):
    """
    Adds 'refine' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser('refine', help='Refine variants.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        '--variant_type',
        type=str,
        required=True,
        choices=VariantTypes.ALL,
        help="Variant type (%s). "
             "If the input VCF file is of structural variants, specify '%s'. "
             "If the input VCF file is of SNVs and INDELs, specify '%s'."
             % (', '.join(f"'{item}'" for item in VariantTypes.ALL),
                VariantTypes.SV,
                VariantTypes.SNV_INDEL)
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
             "This parameter must be specified if --variant_type is '%s'."
             % VariantTypes.SNV_INDEL
    )
    parser_optional.add_argument(
        "--tumor_normal_paired",
        dest="tumor_normal_paired",
        type=bool,
        required=False,
        help="If true, variants were called using tumor and matched normal. "
             "This parameter must be specified if --variant_type is '%s'."
             % VariantTypes.SNV_INDEL
    )
    parser_optional.add_argument(
        "--keep_only_chromosomes",
        dest="keep_only_chromosomes",
        type=str,
        required=False,
        default=[],
        action="append",
        nargs='+',
        help="Chromosomes to keep (e.g. --chromosomes chr1 chr2 chr3). "
             "Chromosomes not specified in this parameter will be removed."
    )
    parser_optional.add_argument(
        "--keep_only_precise_sv",
        dest="keep_only_precise_sv",
        type=bool,
        required=False,
        default=KEEP_ONLY_PRECISE_SV,
        help="Specify as 0 or 1. "
             "If 1 (i.e. true), only retains 'precise' variants (default: 1)."
             " This parameter only applies to structural variants (i.e. --variant_type %s)."
             % VariantTypes.SV
    )
    parser_optional.add_argument(
        "--keep_only_filter_values",
        dest="keep_only_filter_values",
        action="append",
        nargs='+',
        required=False,
        default=KEEP_ONLY_FILTER_VALUES,
        help="VCF 'FILTER' values to keep (default: %s). "
             "Variants that do not have 'FILTER' values specified "
             "in this parameter will be removed."
             % (', '.join(KEEP_ONLY_FILTER_VALUES))
    )
    parser_optional.add_argument(
        "--min_total_coverage",
        dest="min_total_coverage",
        type=int,
        required=False,
        default=MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE,
        help="Minimum total coverage (default: %i)."
             % MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE
    )
    parser_optional.add_argument(
        "--min_variant_reads_count",
        dest="min_variant_reads_count",
        type=int,
        required=False,
        default=MIN_GENOMIC_VARIANT_READS_COUNT,
        help="Minimum variant reads count (default: %i)."
             % MIN_GENOMIC_VARIANT_READS_COUNT
    )
    parser_optional.add_argument(
        "--gapped_regions_tsv_file",
        dest="gapped_regions_tsv_file",
        type=str,
        required=False,
        help="Gapped regions TSV file. "
             "SVs with breakpoints near these gapped regions will be removed. "
             "Expected headers: 'chrom', 'chromStart', 'chromEnd'."
    )
    parser_optional.add_argument(
        "--gapped_regions_padding",
        dest="gapped_regions_padding",
        type=int,
        required=False,
        default=GENOME_GAPPED_REGIONS_PADDING,
        help="Number of bases to pad the gapped regions (default: %i)."
             % GENOME_GAPPED_REGIONS_PADDING
    )
    parser_optional.add_argument(
        "--exclude_sv_tsv_file",
        dest="exclude_sv_tsv_file",
        type=str,
        required=False,
        help="TSV file of SVs to explicitly exclude. "
             "SVs in the input VCF file with breakpoints near the variants "
             "specified in this TSV file will be removed. "
             "Expected headers: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'sv_type'. "
             "This parameter can be used to filter out germline SVs. "
             "Note that this parameter takes into consideration the SV type. "
             "The sv_type must match in order for a putative SV to be filtered out."
    )
    parser_optional.add_argument(
        "--exclude_sv_padding",
        dest="exclude_sv_padding",
        type=int,
        required=False,
        default=EXCLUDE_SV_PADDING,
        help="Number of bases to pad the breakpoints of SVs to exclude "
             "(default: %i)."
             % EXCLUDE_SV_PADDING
    )
    parser_optional.add_argument(
        "--exclude_snv_indel_tsv_file",
        dest="exclude_snv_indel_tsv_file",
        type=str,
        required=False,
        help="TSV file of SNVs and INDELs to explicitly exclude. "
             "SNVs and INDELs in the input VCF file with breakpoints near "
             "the variants specified in this file will be removed. "
             "Expected headers: 'chr_1', 'pos_1', 'chr_2', 'pos_2'. "
             "This parameter can be used to filter out germline SNVs and INDELs. "
             "Keep 'chr_1' the same as 'chr_2' and 'pos_1' the same as 'pos_2' "
             "for SNVs and insertions. Note that this parameter does not take "
             "into consideration the variant type."
    )
    parser_optional.add_argument(
        "--exclude_snv_indel_padding",
        dest="exclude_snv_indel_padding",
        type=int,
        required=False,
        default=EXCLUDE_SNV_INDEL_PADDING,
        help="Number of bases to pad the breakpoints of SNVs and INDELs "
             "to exclude (default: %i)."
             % EXCLUDE_SNV_INDEL_PADDING
    )
    parser.set_defaults(which='refine')
    return sub_parsers


def run_exacto_refine_from_parsed_args(args):
    """
    Run Exacto 'refine' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                variant_type
                vcf_file
                variant_calling_method
                sequencing_platform
                output_tsv_file
                tumor_sample_id
                tumor_normal_paired
                keep_only_chromosomes
                keep_only_precise_sv
                keep_only_filter_values
                min_total_coverage
                min_variant_reads_count
                gapped_regions_tsv_file
                gapped_regions_padding
                exclude_sv_tsv_file
                exclude_sv_padding
                exclude_snv_indel_tsv_file
                exclude_snv_indel_padding
    """
    if args.variant_type == VariantTypes.SV:
        # Step 1. Load VCF file.
        if args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SNIFFLES2:
            df_structural_variants = convert_sniffles2_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.CUTESV:
            df_structural_variants = convert_cutesv_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.SVIM:
            df_structural_variants = convert_svim_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform
            )
        elif args.variant_calling_method == VariantCallingMethods.StructuralVariantCallingMethods.PBSV:
            df_structural_variants = convert_pbsv_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform
            )
        else:
            raise Exception(
                "Invalid value for '--variant_calling_method': %s. "
                "Allowed '--variant_calling_method' values are %s "
                % (args.variant_calling_method,
                   ', '.join(f"'{item}'" for item in VariantCallingMethods.StructuralVariantCallingMethods.ALL))
            )

        # Step 2. Load structural variants to exclude.
        if args.exclude_sv_tsv_file is not None:
            df_structural_variants_to_exclude = pd.read_csv(
                args.exclude_sv_tsv_file,
                sep='\t'
            )
        else:
            df_structural_variants_to_exclude = None

        # Step 3. Load gapped regions.
        if args.gapped_regions_tsv_file is not None:
            df_gapped_regions = pd.read_csv(
                args.gapped_regions_tsv_file,
                sep='\t'
            )
        else:
            df_gapped_regions = None

        # Step 4. Perform refinement.
        df_structural_variants = run_exacto_refine_genomic_structural_variants(
            df_structural_variants=df_structural_variants,
            df_structural_variants_to_exclude=df_structural_variants_to_exclude,
            df_gapped_regions=df_gapped_regions,
            variant_calling_method=args.variant_calling_method,
            keep_only_precise_sv=args.keep_only_precise_sv,
            keep_only_chromosomes=args.keep_only_chromosomes,
            keep_only_filter_values=args.keep_only_filter_values,
            min_total_coverage=args.min_total_coverage,
            min_variant_reads_count=args.min_variant_reads_count,
            gapped_regions_padding=args.gapped_regions_padding,
            exclude_variants_padding=args.exclude_sv_padding
        )

        # Step 5. Write refined structural variants to a TSV file.
        df_structural_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
    elif args.variant_type == VariantTypes.SNV_INDEL:
        # Step 1. Load VCF file.
        if args.variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.GATK4_MUTECT2:
            df_variants = convert_gatk4_mutect2_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform,
                tumor_sample_id=args.tumor_sample_id
            )
        elif args.variant_calling_method == VariantCallingMethods.SmallVariantCallingMethods.DEEPVARIANT:
            df_variants = convert_deepvariant_vcf_to_dataframe(
                vcf_file=args.vcf_file,
                sequencing_platform=args.sequencing_platform
            )
        else:
            raise Exception(
                "Invalid value for '--variant_calling_method': %s. "
                "Allowed '--variant_calling_method' values are %s "
                % (args.variant_calling_method,
                   ', '.join(f"'{item}'" for item in VariantCallingMethods.SmallVariantCallingMethods.ALL))
            )
        # Step 2. Load gapped regions.
        if args.gapped_regions_tsv_file is not None:
            df_gapped_regions = pd.read_csv(
                args.gapped_regions_tsv_file,
                sep='\t'
            )
        else:
            df_gapped_regions = None

        # Step 3. Perform refinement.
        df_variants = run_exacto_refine_genomic_small_variants(
            df_variants=df_variants,
            df_gapped_regions=df_gapped_regions,
            variant_calling_method=args.variant_calling_method,
            tumor_normal_paired=args.tumor_normal_paired,
            keep_only_chromosomes=args.keep_only_chromosomes,
            keep_only_filter_values=args.keep_only_filter_values,
            min_total_coverage=args.min_total_coverage,
            min_variant_reads_count=args.min_variant_reads_count,
            gapped_regions_padding=args.gapped_regions_padding
        )

        # Step 4. Write refined variants to a TSV file.
        df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)
    else:
        raise Exception(
            "Invalid value for '--variant_type': %s. "
            "Allowed '--variant_type' values are %s "
            % (args.variant_type,
               ', '.join(f"'{item}'" for item in VariantTypes.ALL))
        )
