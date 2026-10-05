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
and run Exacto 'stitch-reference-transcripts' command.
"""


import argparse
import os

from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_stitch_reference_transcripts_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'stitch-reference-transcripts' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser(
        'stitch-reference-transcripts',
        help="Stitch reference transcript sequence onto the degraded ends of aligned RNA sequences. "
             "The stitched sequences no longer match the read coordinates of any "
             "call-rna-transcript-vars output computed on the same BAM file: realign the stitched "
             "FASTA and re-run call-rna-transcript-vars before translate-transcripts.")
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file of the aligned RNA sequences to stitch (e.g. the realigned consensus "
             "sequences of exacto determine-rna-consensus). Every mapped read is one transcript, "
             "identified by its read name."
    )
    parser_required.add_argument(
        "--reference-genome-fasta-file",
        dest="reference_genome_fasta_file",
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-file",
        dest="reference_gene_annotation_file",
        type=str,
        required=True,
        help="Reference gene annotation file."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-source",
        dest="reference_gene_annotation_source",
        type=str,
        required=True,
        help="Reference gene annotation source (choices: %s)." %
             ','.join([str(GeneAnnotationSource.GENCODE)])
    )
    parser_required.add_argument(
        "--reference-gene-annotation-assembly",
        dest="reference_gene_annotation_assembly",
        type=str,
        required=True,
        help="Reference gene annotation assembly (e.g. 'hg38')."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-version",
        dest="reference_gene_annotation_version",
        type=str,
        required=True,
        help="Reference gene annotation version (e.g. 'v41')."
    )
    parser_required.add_argument(
        "--output-dir",
        dest="output_dir",
        type=str,
        required=True,
        help="Output directory."
    )
    parser_required.add_argument(
        "--output-prefix",
        dest="output_prefix",
        type=str,
        required=True,
        help="Output prefix."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality. A read whose best-mapped record falls below this is left "
             "out of the output (default: %i)." % STITCH_REFERENCE_TRANSCRIPTS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-num-splice-junction-matches",
        dest="min_num_splice_junction_matches",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_MIN_NUM_SPLICE_JUNCTION_MATCHES,
        required=False,
        help="A reference transcript match with fewer splice junction matches than this is "
             "treated as no match. A single-exon reference transcript is exempt "
             "(default: %i)." % STITCH_REFERENCE_TRANSCRIPTS_MIN_NUM_SPLICE_JUNCTION_MATCHES
    )
    parser_optional.add_argument(
        "--pas-search-size",
        dest="pas_search_size",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_PAS_SEARCH_SIZE,
        required=False,
        help="Number of read bases, ending at an intronic 3' terminal aligned base, searched for a "
             "polyadenylation signal hexamer (default: %i)."
             % STITCH_REFERENCE_TRANSCRIPTS_PAS_SEARCH_SIZE
    )
    parser_optional.add_argument(
        "--pas-hexamers",
        dest="pas_hexamers",
        type=str,
        nargs='+',
        default=STITCH_REFERENCE_TRANSCRIPTS_PAS_HEXAMERS,
        required=False,
        help="Polyadenylation signal hexamers, in transcript orientation. Pass multiple as: "
             "--pas-hexamers AATAAA ATTAAA (default: %s)."
             % ' '.join(STITCH_REFERENCE_TRANSCRIPTS_PAS_HEXAMERS)
    )
    parser_optional.add_argument(
        "--pas-start-offset-range",
        dest="pas_start_offset_range",
        type=int,
        nargs=2,
        metavar=('MIN', 'MAX'),
        default=STITCH_REFERENCE_TRANSCRIPTS_PAS_START_OFFSET_RANGE,
        required=False,
        help="Permitted distance in bases, inclusive, from a polyadenylation signal hexamer's first "
             "base to the end of the searched read sequence. Pass as: --pas-start-offset-range 10 40 "
             "(default: %i %i)." % tuple(STITCH_REFERENCE_TRANSCRIPTS_PAS_START_OFFSET_RANGE)
    )
    parser_optional.add_argument(
        "--polya-tail-min-adenosine-fraction",
        dest="polya_tail_min_adenosine_fraction",
        type=float,
        default=STITCH_REFERENCE_TRANSCRIPTS_POLYA_TAIL_MIN_ADENOSINE_FRACTION,
        required=False,
        help="Read bases past the 3' terminal aligned base are accepted as a polyA tail when at "
             "least this fraction of them is A; otherwise they may be the remnant of a fusion or "
             "aberrant splicing, and the 3' end is left alone (default: %s)."
             % STITCH_REFERENCE_TRANSCRIPTS_POLYA_TAIL_MIN_ADENOSINE_FRACTION
    )
    parser_optional.add_argument(
        "--polya-window-size",
        dest="polya_window_size",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_POLYA_WINDOW_SIZE,
        required=False,
        help="Number of genomic bases on each side of an intronic 3' terminal aligned base scanned "
             "for an A-tract (default: %i)." % STITCH_REFERENCE_TRANSCRIPTS_POLYA_WINDOW_SIZE
    )
    parser_optional.add_argument(
        "--polya-min-consecutive-adenosine",
        dest="polya_min_consecutive_adenosine",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_POLYA_MIN_CONSECUTIVE_ADENOSINE,
        required=False,
        help="An A-run of at least this length in that genomic window calls internal priming. An "
             "intronic 3' end with a polyadenylation signal and no such A-run is kept as a genuine "
             "intronic polyadenylation site and is not stitched (default: %i)."
             % STITCH_REFERENCE_TRANSCRIPTS_POLYA_MIN_CONSECUTIVE_ADENOSINE
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=STITCH_REFERENCE_TRANSCRIPTS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % STITCH_REFERENCE_TRANSCRIPTS_NUM_THREADS
    )
    parser.set_defaults(which='stitch-reference-transcripts')
    return sub_parsers


def run_cli_stitch_reference_transcripts_from_parsed_args(args) -> None:
    """
    Run Exacto 'stitch-reference-transcripts' command using parameters from parsed arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    stitch_reference_transcripts(
        bam_file=args.bam_file,
        reference_genome_fasta_file=args.reference_genome_fasta_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_source=GeneAnnotationSource(args.reference_gene_annotation_source),
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        min_mapping_quality=args.min_mapping_quality,
        min_num_splice_junction_matches=args.min_num_splice_junction_matches,
        pas_search_size=args.pas_search_size,
        pas_hexamers=args.pas_hexamers,
        pas_start_offset_range=tuple(args.pas_start_offset_range),
        polya_tail_min_adenosine_fraction=args.polya_tail_min_adenosine_fraction,
        polya_window_size=args.polya_window_size,
        polya_min_consecutive_adenosine=args.polya_min_consecutive_adenosine,
        num_threads=args.num_threads
    )
