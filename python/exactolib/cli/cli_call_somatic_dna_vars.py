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
and run Exacto 'call-somatic-dna-vars' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_call_somatic_dna_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'call-somatic-dna-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('call-somatic-dna-vars', help='Call somatic DNA variants in a long-read WGS BAM file.')
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
        "--control-bam-files",
        dest="control_bam_files",
        type=str,
        nargs='+',
        default=[],
        required=True,
        help="Input control BAM file(s) (e.g. --control-bam-files BAM_FILE_1 BAM_FILE_2)."
    )
    parser_required.add_argument(
        "--control-bam-bai-files",
        dest="control_bam_bai_files",
        type=str,
        nargs='+',
        default=[],
        required=True,
        help="Input control BAM.BAI file(s) (e.g. --control-bam-bai-files BAM_BAI_FILE_1 BAM_BAI_FILE_2)."
    )
    parser_required.add_argument(
        "--fasta-file",
        dest="fasta_file",
        type=str,
        required=True,
        help="Input reference genome FASTA file."
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
        "--preset",
        dest="preset",
        type=str,
        choices=list(CALL_DNA_VARS_PRESETS.keys()),
        default=None,
        required=False,
        help="Sequencing-platform preset that fills platform-typical defaults "
             "for --expected-sequencing-error, and --expected-slippage-probability. Choices: 'pb' (PacBio HiFi) or "
             "'ont' (Oxford Nanopore). Any explicit parameter wins over the preset."
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_SOMATIC_DNA_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        '--regions',
        dest='regions',
        type=str,
        nargs='+',
        default=[],
        required=False,
        help='Genomic regions in which to identify variants (e.g. --regions chr1 chr2 or --regions chr1:1-1000000 chr2:1-1000000). '
             'If unspecified, Exacto identifies variants in all contigs '
             'found in the BAM file (--bam-file BAM_FILE).'
    )
    parser_optional.add_argument(
        "--min-reads",
        dest="min_reads",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_READS,
        required=False,
        help="Minimum number of supporting reads (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_READS
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-base-quality",
        dest="min_base_quality",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_BASE_QUALITY,
        required=False,
        help="Minimum base quality (default: %f)."
             % CALL_SOMATIC_DNA_VARS_MIN_BASE_QUALITY
    )
    parser_optional.add_argument(
        "--min-total-depth",
        dest="min_total_depth",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_TOTAL_DEPTH,
        required=False,
        help="Minimum total depth (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_TOTAL_DEPTH
    )
    parser_optional.add_argument(
        "--min-alt-allele-fraction",
        dest="min_alt_allele_fraction",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_MIN_ALT_ALLELE_FRACTION,
        required=False,
        help="Minimum alternate allele fraction (default: %f)."
             % CALL_SOMATIC_DNA_VARS_MIN_ALT_ALLELE_FRACTION
    )
    parser_optional.add_argument(
        "--min-size-proportion",
        dest="min_size_proportion",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_MIN_SIZE_PROPORTION,
        required=False,
        help="Minimum size proportion between two variants (default: %f). "
             "Size proportion = smaller variant size / longer variant size."
             % CALL_SOMATIC_DNA_VARS_MIN_SIZE_PROPORTION
    )
    parser_optional.add_argument(
        "--max-ins-norm-edit-distance",
        dest="max_ins_norm_edit_distance",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE,
        required=False,
        help="Maximum insertion normalized edit (Levenshtein) distance (default: %f). "
             "Normalized edit distance = edit distance / longer insertion size."
             % CALL_SOMATIC_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE
    )
    parser_optional.add_argument(
        "--max-intrachromosomal-distance",
        dest="max_intrachromosomal_distance",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE,
        required=False,
        help="Maximum distance for clustering intrachromomsomal variants (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE
    )
    parser_optional.add_argument(
        "--max-intrachromosomal-distance-tau",
        dest="max_intrachromosomal_distance_tau",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU,
        required=False,
        help="Maximum distance tau for clustering intrachromomsomal variants (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU
    )
    parser_optional.add_argument(
        "--max-interchromosomal-distance",
        dest="max_interchromosomal_distance",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE,
        required=False,
        help="Maximum distance for clustering intrachromomsomal variants (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE
    )
    parser_optional.add_argument(
        "--max-slippage-repeat-length",
        dest="max_slippage_repeat_length",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH,
        required=False,
        help="Maximum slippage repeat length (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MAX_SLIPPAGE_REPEAT_LENGTH
    )
    parser_optional.add_argument(
        "--apply-infinite-sites-assumption",
        dest="apply_infinite_sites_assumption",
        type=str2bool,
        default=CALL_SOMATIC_DNA_VARS_INFINITE_SITES_ASSUMPTION,
        required=False,
        help="If 'yes', apply infinite sites assumption to the variant calling. That is, "
             "if a variant in the BAM file shares breakpoint with any of the variant in"
             "any of the control BAM files, filter it out (default: %s)."
             % CALL_SOMATIC_DNA_VARS_INFINITE_SITES_ASSUMPTION
    )
    parser_optional.add_argument(
        "--chunk-size",
        dest="chunk_size",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_CHUNK_SIZE,
        required=False,
        help="Chunk size for variant calling (default: %i)." % CALL_SOMATIC_DNA_VARS_CHUNK_SIZE
    )
    parser_optional.add_argument(
        "--max-records",
        dest="max_records",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_RECORDS,
        required=False,
        help="Maximum number of records. Read names having more than this value will be excluded (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MAX_RECORDS
    )
    parser_optional.add_argument(
        "--expected-variant-allele-fraction",
        dest="expected_variant_allele_fraction",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION,
        required=False,
        help="Expected variant allele fraction (default: %f)."
             % CALL_SOMATIC_DNA_VARS_EXPECTED_VARIANT_ALLELE_FRACTION
    )
    parser_optional.add_argument(
        "--expected-mutation-rate",
        dest="expected_mutation_rate",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_EXPECTED_MUTATION_RATE,
        required=False,
        help="Expected mutation rate (default: %f)."
             % CALL_SOMATIC_DNA_VARS_EXPECTED_MUTATION_RATE
    )
    parser_optional.add_argument(
        "--expected-sequencing-error",
        dest="expected_sequencing_error",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
        required=False,
        help="Expected sequencing error rate (default: %f)."
             % CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR
    )
    parser_optional.add_argument(
        "--expected-slippage-probability",
        dest="expected_slippage_probability",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
        required=False,
        help="Expected slippage probability (default: %f). Suggested: 2x the expected sequencing error rate."
             % CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY
    )
    parser_optional.add_argument(
        "--max-f1-fraction",
        dest="max_f1_fraction",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_MAX_F1_FRACTION,
        required=False,
        help="Maximum F1 value fraction (default: %f)."
             % CALL_SOMATIC_DNA_VARS_MAX_F1_FRACTION
    )
    parser_optional.add_argument(
        "--max-fpr",
        dest="max_fpr",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_MAX_FPR,
        required=False,
        help="Maximum false positive rate (default: %f)."
             % CALL_SOMATIC_DNA_VARS_MAX_FPR
    )
    parser_optional.add_argument(
        "--temp-dir",
        dest="temp_dir",
        type=str,
        default="",
        required=False,
        help="Temp directory (default: TMPDIR)."
    )
    parser.set_defaults(which='call-somatic-dna-vars')
    return sub_parsers


def run_cli_call_somatic_dna_vars_from_parsed_args(args):
    """
    Run Exacto 'call-somatic-dna-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    bam_file
                    bam_bai_file
                    mode
                    output_tsv_file
                    control_bam_file
                    control_bam_bai_file
                    num_threads
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
    if len(args.regions) == 0:
        regions = []
    else:
        chromosome_lengths = get_chromosome_lengths(bam_file=args.bam_file)
        regions = []
        for region in args.regions:
            if ':' in region:
                chromosome = region.split(':')[0]
                start = int(region.split(':')[1].split('-')[0])
                end = int(region.split(':')[1].split('-')[1])
            else:
                chromosome = region
                start = 1
                end = chromosome_lengths[chromosome]
            regions.append((chromosome, start, end))

    identify_somatic_dna_variants(
        bam_file=args.bam_file,
        bam_bai_file=args.bam_bai_file,
        control_bam_files=args.control_bam_files,
        control_bam_bai_files=args.control_bam_bai_files,
        fasta_file=args.fasta_file,
        output_tsv_file=args.output_tsv_file,
        regions=regions,
        min_reads=args.min_reads,
        min_mapping_quality=args.min_mapping_quality,
        min_base_quality=args.min_base_quality,
        min_total_depth=args.min_total_depth,
        min_alt_allele_fraction=args.min_alt_allele_fraction,
        min_size_proportion=args.min_size_proportion,
        max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
        max_intrachromosomal_distance=args.max_intrachromosomal_distance,
        max_intrachromosomal_distance_tau=args.max_intrachromosomal_distance_tau,
        max_interchromosomal_distance=args.max_interchromosomal_distance,
        max_slippage_repeat_length=args.max_slippage_repeat_length,
        num_threads=args.num_threads,
        chunk_size=args.chunk_size,
        max_records=args.max_records,
        expected_variant_allele_fraction=args.expected_variant_allele_fraction,
        expected_mutation_rate=args.expected_mutation_rate,
        expected_sequencing_error=args.expected_sequencing_error,
        expected_slippage_probability=args.expected_slippage_probability,
        max_f1_fraction=args.max_f1_fraction,
        apply_infinite_sites_assumption=args.apply_infinite_sites_assumption,
        temp_dir=args.temp_dir
    )
