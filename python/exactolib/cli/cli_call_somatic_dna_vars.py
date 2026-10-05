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
        "--bai-file",
        dest="bai_file",
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
        "--control-bai-files",
        dest="control_bai_files",
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
             "for --expected-sequencing-error, --expected-slippage-probability, --min-reads, and "
             "--min-total-depth. Choices: 'pb' (PacBio HiFi) or "
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
        default=None,
        required=False,
        help="Minimum number of supporting reads (default: preset-dependent, pb 4 / ont 5; "
             "%i with no preset)."
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
        help="Minimum base quality (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_BASE_QUALITY
    )
    parser_optional.add_argument(
        "--min-total-depth",
        dest="min_total_depth",
        type=int,
        default=None,
        required=False,
        help="Minimum total depth (default: preset-dependent, pb 4 / ont 5; %i with no preset)."
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
        default=None,
        required=False,
        help="Expected sequencing error rate (default: %f with no preset, %f with --preset pb, %f with --preset ont)."
             % (CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR,
                CALL_DNA_VARS_PRESETS["pb"]["expected_sequencing_error"],
                CALL_DNA_VARS_PRESETS["ont"]["expected_sequencing_error"])
    )
    parser_optional.add_argument(
        "--expected-slippage-probability",
        dest="expected_slippage_probability",
        type=float,
        default=None,
        required=False,
        help="Expected slippage probability (default: %f with no preset, %f with --preset pb, %f with --preset ont). "
             "Suggested: 2x the expected sequencing error rate."
             % (CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY,
                CALL_DNA_VARS_PRESETS["pb"]["expected_slippage_probability"],
                CALL_DNA_VARS_PRESETS["ont"]["expected_slippage_probability"])
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
    parser_optional.add_argument(
        "--min-terminal-softclip-length",
        dest="min_terminal_softclip_length",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        required=False,
        help="Terminal soft-clip runs shorter than this many bases emit no variant record; 1-2 bp "
             "terminal clips are aligner anchor jitter, not junction or insertion evidence. Set 0 to "
             "disable. Currently has no effect: the somatic pipeline always uses 0 and leaves clip "
             "evidence to breakend clustering (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH
    )
    parser_optional.add_argument(
        "--bkpt-rescue",
        dest="bkpt_rescue",
        type=str2bool,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE,
        required=False,
        help="Whether to rescue breakpoints from insertions by re-aligning long insertion sequences "
             "(default: %s)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-min-ins-len",
        dest="bkpt_rescue_min_ins_len",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MIN_INS_LEN,
        required=False,
        help="Minimum insertion length to be considered for re-alignment (to rescue breakpoints) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MIN_INS_LEN
    )
    parser_optional.add_argument(
        "--bkpt-rescue-max-ins-len",
        dest="bkpt_rescue_max_ins_len",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MAX_INS_LEN,
        required=False,
        help="Maximum insertion length to be considered for re-alignment (to rescue breakpoints) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_MAX_INS_LEN
    )
    parser_optional.add_argument(
        "--bkpt-rescue-search-distance",
        dest="bkpt_rescue_search_distance",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_SEARCH_DISTANCE,
        required=False,
        help="Search distance in bases around an insertion for re-alignment (to rescue breakpoints) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_SEARCH_DISTANCE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-open-score",
        dest="bkpt_rescue_realignment_gap_open_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        required=False,
        help="Gap open score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-extend-score",
        dest="bkpt_rescue_realignment_gap_extend_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        required=False,
        help="Gap extend score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-k",
        dest="bkpt_rescue_realignment_k",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_K,
        required=False,
        help="K-mer size for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_K
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-band-width",
        dest="bkpt_rescue_realignment_band_width",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        required=False,
        help="Band width for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-score-fraction",
        dest="bkpt_rescue_realignment_min_score_fraction",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        required=False,
        help="Minimum score fraction for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-query-coverage",
        dest="bkpt_rescue_realignment_min_query_coverage",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        required=False,
        help="Minimum query coverage for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-span-proportion",
        dest="bkpt_rescue_realignment_min_span_proportion",
        type=float,
        default=CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SPAN_PROPORTION,
        required=False,
        help="Minimum span proportion for insertion re-alignment (to rescue breakpoints): the shorter "
             "of the insertion length and the reference span its re-aligned sequence covers, over the "
             "longer (default: %f)."
             % CALL_SOMATIC_DNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SPAN_PROPORTION
    )
    parser_optional.add_argument(
        "--poa-match-score",
        dest="poa_match_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_POA_MATCH_SCORE,
        required=False,
        help="Partial order alignment match score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_POA_MATCH_SCORE
    )
    parser_optional.add_argument(
        "--poa-mismatch-score",
        dest="poa_mismatch_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_POA_MISMATCH_SCORE,
        required=False,
        help="Partial order alignment mismatch score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_POA_MISMATCH_SCORE
    )
    parser_optional.add_argument(
        "--poa-gap-open-score",
        dest="poa_gap_open_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_POA_GAP_OPEN_SCORE,
        required=False,
        help="Partial order alignment gap open score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_POA_GAP_OPEN_SCORE
    )
    parser_optional.add_argument(
        "--poa-gap-extend-score",
        dest="poa_gap_extend_score",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_POA_GAP_EXTEND_SCORE,
        required=False,
        help="Partial order alignment gap extend score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CALL_SOMATIC_DNA_VARS_POA_GAP_EXTEND_SCORE
    )
    parser_optional.add_argument(
        "--read-depth-max-merge-distance",
        dest="read_depth_max_merge_distance",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_READ_DEPTH_MAX_MERGE_DISTANCE,
        required=False,
        help="Positions less than this many bases apart have their read depths counted from one query "
             "of the BAM file. Only the number of queries depends on it, not the depths (default: %i)."
             % CALL_SOMATIC_DNA_VARS_READ_DEPTH_MAX_MERGE_DISTANCE
    )
    parser_optional.add_argument(
        "--min-homopolymer-len",
        dest="min_homopolymer_len",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_HOMOPOLYMER_LEN,
        required=False,
        help="Minimum homopolymer length for a variant to be treated as homopolymer-context (default: "
             "%i)."
             % CALL_SOMATIC_DNA_VARS_MIN_HOMOPOLYMER_LEN
    )
    parser_optional.add_argument(
        "--min-dinucleotide-context-len",
        dest="min_dinucleotide_context_len",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MIN_DINUCLEOTIDE_CONTEXT_LEN,
        required=False,
        help="Minimum dinucleotide repeat context length for a variant to be treated as "
             "dinucleotide-context (default: %i)."
             % CALL_SOMATIC_DNA_VARS_MIN_DINUCLEOTIDE_CONTEXT_LEN
    )
    parser_optional.add_argument(
        "--max-control-reads",
        dest="max_control_reads",
        type=int,
        default=CALL_SOMATIC_DNA_VARS_MAX_CONTROL_READS,
        required=False,
        help="Maximum number of control reads a variant may have and still be called somatic (default: "
             "%i)."
             % CALL_SOMATIC_DNA_VARS_MAX_CONTROL_READS
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
                    temp_dir
    """
    # Resolve preset-affected parameters: explicit flag > preset > built-in default.
    preset_values: dict[str, int | float] = CALL_DNA_VARS_PRESETS.get(args.preset, {})
    expected_sequencing_error = (
        args.expected_sequencing_error
        if args.expected_sequencing_error is not None
        else preset_values.get("expected_sequencing_error", CALL_SOMATIC_DNA_VARS_EXPECTED_SEQUENCING_ERROR)
    )
    expected_slippage_probability = (
        args.expected_slippage_probability
        if args.expected_slippage_probability is not None
        else preset_values.get("expected_slippage_probability", CALL_SOMATIC_DNA_VARS_EXPECTED_SLIPPAGE_PROBABILITY)
    )
    min_reads = (
        args.min_reads
        if args.min_reads is not None
        else preset_values.get("min_reads", CALL_SOMATIC_DNA_VARS_MIN_READS)
    )
    min_total_depth = (
        args.min_total_depth
        if args.min_total_depth is not None
        else preset_values.get("min_total_depth", CALL_SOMATIC_DNA_VARS_MIN_TOTAL_DEPTH)
    )

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
        bam_bai_file=args.bai_file,
        control_bam_files=args.control_bam_files,
        control_bam_bai_files=args.control_bai_files,
        fasta_file=args.fasta_file,
        output_tsv_file=args.output_tsv_file,
        regions=regions,
        min_reads=min_reads,
        min_mapping_quality=args.min_mapping_quality,
        min_base_quality=args.min_base_quality,
        min_total_depth=min_total_depth,
        min_alt_allele_fraction=args.min_alt_allele_fraction,
        min_size_proportion=args.min_size_proportion,
        max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
        max_intrachromosomal_distance=args.max_intrachromosomal_distance,
        max_slippage_repeat_length=args.max_slippage_repeat_length,
        min_terminal_softclip_length=args.min_terminal_softclip_length,
        bkpt_rescue=args.bkpt_rescue,
        bkpt_rescue_min_ins_len=args.bkpt_rescue_min_ins_len,
        bkpt_rescue_max_ins_len=args.bkpt_rescue_max_ins_len,
        bkpt_rescue_search_distance=args.bkpt_rescue_search_distance,
        bkpt_rescue_realignment_gap_open_score=args.bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=args.bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=args.bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=args.bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=args.bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=args.bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_span_proportion=args.bkpt_rescue_realignment_min_span_proportion,
        poa_match_score=args.poa_match_score,
        poa_mismatch_score=args.poa_mismatch_score,
        poa_gap_open_score=args.poa_gap_open_score,
        poa_gap_extend_score=args.poa_gap_extend_score,
        read_depth_max_merge_distance=args.read_depth_max_merge_distance,
        min_homopolymer_len=args.min_homopolymer_len,
        min_dinucleotide_context_len=args.min_dinucleotide_context_len,
        max_control_reads=args.max_control_reads,
        num_threads=args.num_threads,
        chunk_size=args.chunk_size,
        max_records=args.max_records,
        expected_variant_allele_fraction=args.expected_variant_allele_fraction,
        expected_mutation_rate=args.expected_mutation_rate,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_fpr=args.max_fpr,
        apply_infinite_sites_assumption=args.apply_infinite_sites_assumption,
        temp_dir=args.temp_dir
    )
