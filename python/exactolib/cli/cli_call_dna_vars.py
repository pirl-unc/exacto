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
and run Exacto 'call-dna-vars' command.
"""


import argparse
from ..constants import *
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_call_dna_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'call-dna-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('call-dna-vars', help='Call DNA variants in a long-read WGS BAM file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file."
    )
    parser_required.add_argument(
        "--bam-bai-file",
        dest="bam_bai_file",
        type=str,
        required=True,
        help="Input BAM.BAI file."
    )
    parser_required.add_argument(
        "--mode",
        dest="mode",
        type=str,
        choices=[DnaVariantCallingModes.ALL, DnaVariantCallingModes.CASE_SPECIFIC],
        required=True,
        help="DNA variant calling mode (either 'case-specific' or 'all'). "
             "If --mode case-specific, then at least one control BAM file needs to be supplied."
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
        "--control-bam-files",
        dest="control_bam_files",
        type=str,
        nargs='+',
        default=[],
        required=False,
        help="Input control BAM file(s) (e.g. --control-bam-files BAM_FILE_1 BAM_FILE_2)."
    )
    parser_optional.add_argument(
        "--control-bam-bai-files",
        dest="control_bam_bai_files",
        type=str,
        nargs='+',
        default=[],
        required=False,
        help="Input control BAM.BAI file(s) (e.g. --control-bam-bai-files BAM_BAI_FILE_1 BAM_BAI_FILE_2)."
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_DNA_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_DNA_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        default=CALL_DNA_VARS_GZIP,
        required=False,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % CALL_DNA_VARS_GZIP
    )
    parser_optional.add_argument(
        '--chromosomes',
        dest='chromosomes',
        type=str,
        nargs='+',
        default=[],
        required=False,
        help='Chromosomes in which to identify variants (e.g. --chromosomes chr1 chr2). '
             'If unspecified, Exacto identifies variants in all chromosomes '
             'found in the BAM file (--bam-file BAM_FILE).'
    )
    parser_optional.add_argument(
        "--min-reads",
        dest="min_reads",
        type=int,
        default=CALL_DNA_VARS_MIN_READS,
        required=False,
        help="Minimum number of supporting reads (default: %i)."
             % CALL_DNA_VARS_MIN_READS
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=CALL_DNA_VARS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % CALL_DNA_VARS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-average-base-quality",
        dest="min_average_base_quality",
        type=float,
        default=CALL_DNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        required=False,
        help="Minimum average base quality (default: %f)."
             % CALL_DNA_VARS_MIN_AVERAGE_BASE_QUALITY
    )
    parser_optional.add_argument(
        "--min-size-proportion",
        dest="min_size_proportion",
        type=float,
        default=CALL_DNA_VARS_MIN_SIZE_PROPORTION,
        required=False,
        help="Minimum size proportion between two variants (default: %f). "
             "Size proportion = smaller variant size / longer variant size."
             % CALL_DNA_VARS_MIN_SIZE_PROPORTION
    )
    parser_optional.add_argument(
        "--max-ins-norm-edit-distance",
        dest="max_ins_norm_edit_distance",
        type=float,
        default=CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        required=False,
        help="Maximum insertion normalized edit (Levenshtein) distance (default: %f). "
             "Normalized edit distance = edit distance / longer insertion size."
             % CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE
    )
    parser_optional.add_argument(
        "--max-intrachromosomal-distance",
        dest="max_intrachromosomal_distance",
        type=int,
        default=CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        required=False,
        help="Maximum distance for clustering intrachromomsomal variants (default: %i)."
             % CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE
    )
    parser_optional.add_argument(
        "--max-intrachromosomal-distance-tau",
        dest="max_intrachromosomal_distance_tau",
        type=int,
        default=CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        required=False,
        help="Maximum distance tau for clustering intrachromomsomal variants (default: %i)."
             % CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU
    )
    parser_optional.add_argument(
        "--max-interchromosomal-distance",
        dest="max_interchromosomal_distance",
        type=int,
        default=CALL_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        required=False,
        help="Maximum distance for clustering intrachromomsomal variants (default: %i)."
             % CALL_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE
    )
    parser_optional.add_argument(
        "--apply-infinite-sites-assumption",
        dest="apply_infinite_sites_assumption",
        type=str2bool,
        default=CALL_DNA_VARS_INFINITE_SITES_ASSUMPTION,
        required=False,
        help="If 'yes', apply infinite sites assumption to the variant calling. That is, "
             "if a variant in the case BAM file shares breakpoint with any of the variant in"
             "any of the control BAM files, filter it out (default: %s)."
             % CALL_DNA_VARS_INFINITE_SITES_ASSUMPTION
    )
    parser_optional.add_argument(
        "--temp-dir",
        dest="temp_dir",
        type=str,
        default="",
        required=False,
        help="Temp directory (default: TMPDIR)."
    )
    parser.set_defaults(which='call-dna-vars')
    return sub_parsers


def run_cli_call_dna_vars_from_parsed_args(args):
    """
    Run Exacto 'call-dna-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    bam_file
                    bam_bai_file
                    mode
                    output_tsv_file
                    control_bam_file
                    control_bam_bai_file
                    num_threads
                    gzip
                    chromosomes
                    min_reads
                    min_mapping_quality
                    min_average_base_quality
                    min_size_proportion
                    max_ins_norm_edit_distance
                    max_intrachromosomal_distance
                    max_intrachromosomal_distance_tau
                    max_interchromosomal_distance
                    temp_dir
    """
    if len(args.chromosomes) == 0 :
        args.chromosomes = get_chromosomes(bam_file=args.bam_file)

    if args.gzip:
        if args.output_tsv_file.endswith('.gz'):
            output_tsv_file = args.output_tsv_file
        else:
            output_tsv_file = args.output_tsv_file + '.gz'
    else:
        output_tsv_file = args.output_tsv_file

    if args.mode == DnaVariantCallingModes.ALL:
        identify_dna_variants(
            bam_file=args.bam_file,
            bam_bai_file=args.bam_bai_file,
            output_tsv_file=output_tsv_file,
            gzip=args.gzip,
            min_reads=args.min_reads,
            min_mapping_quality=args.min_mapping_quality,
            min_average_base_quality=args.min_average_base_quality,
            min_size_proportion=args.min_size_proportion,
            max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
            max_intrachromosomal_distance=args.max_intrachromosomal_distance,
            max_intrachromosomal_distance_tau=args.max_intrachromosomal_distance_tau,
            max_interchromosomal_distance=args.max_interchromosomal_distance,
            num_threads=args.num_threads,
            chromosomes=args.chromosomes,
            temp_dir=args.temp_dir
        )
    elif args.mode == DnaVariantCallingModes.CASE_SPECIFIC:
        assert(len(args.control_bam_files) > 0)
        assert (len(args.control_bam_files) == len(args.control_bam_bai_files))
        identify_case_specific_dna_variants(
            case_bam_file=args.bam_file,
            case_bam_bai_file=args.bam_bai_file,
            control_bam_files=args.control_bam_files,
            control_bam_bai_files=args.control_bam_bai_files,
            output_tsv_file=output_tsv_file,
            gzip=args.gzip,
            chromosomes=args.chromosomes,
            min_reads=args.min_reads,
            min_mapping_quality=args.min_mapping_quality,
            min_average_base_quality=args.min_average_base_quality,
            min_size_proportion=args.min_size_proportion,
            max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
            max_intrachromosomal_distance=args.max_intrachromosomal_distance,
            max_intrachromosomal_distance_tau=args.max_intrachromosomal_distance_tau,
            max_interchromosomal_distance=args.max_interchromosomal_distance,
            apply_infinite_sites_assumption=args.apply_infinite_sites_assumption,
            num_threads=args.num_threads,
            temp_dir=args.temp_dir
        )
    else:
        raise Exception('Unsupported mode: %s' % args.mode)

