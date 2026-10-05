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
and run Exacto 'translate-transcripts' command.
"""


import argparse
import os
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_translate_transcripts_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'translate-transcripts' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('translate-transcripts', help='Translate transcripts to proteoforms.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--assembled-transcript-model-alignments-tsv-file",
        dest="assembled_transcript_model_alignments_tsv_file",
        type=str,
        required=True,
        help="Input assembled transcript model alignments TSV file (output from exacto call-rna-transcript-vars)."
    )
    parser_required.add_argument(
        "--assembled-transcript-variants-tsv-file",
        dest="assembled_transcript_variants_tsv_file",
        type=str,
        required=True,
        help="Input assembled transcript variants TSV file (output from exacto call-rna-transcript-vars)."
    )
    parser_required.add_argument(
        "--dna-variants-tsv-file",
        dest="dna_variants_tsv_file",
        type=str,
        required=True,
        help="Input DNA variants TSV file (output from exacto call-somatic-dna-vars)."
    )
    parser_required.add_argument(
        "--integrated-variants-tsv-file",
        dest="integrated_variants_tsv_file",
        type=str,
        required=True,
        help="Input integrated variants TSV file (output from exacto integrate-vars)."
    )
    parser_required.add_argument(
        "--strategy",
        dest="strategy",
        type=str,
        choices=[str(TranslationStrategy.LONGEST_ORF), str(TranslationStrategy.ALL_ORFS)],
        required=True,
        help="Translation strategy."
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
        default="",
        help="Output prefix."
    )

    # Transcript sequences and their supporting reads: exactly one source.
    parser_sequences = parser.add_argument_group(
        'transcript sequence arguments (exactly one required)'
    )
    parser_sequences_group = parser_sequences.add_mutually_exclusive_group(required=True)
    parser_sequences_group.add_argument(
        "--rna-consensus-tsv-file",
        dest="rna_consensus_tsv_file",
        type=str,
        default="",
        help="Input RNA consensus TSV file (output from exacto determine-rna-consensus; "
             "columns: cluster_id, consensus_sequence, num_reads, read_names). Transcripts are "
             "named by cluster ID, matching the FASTA that command writes."
    )
    parser_sequences_group.add_argument(
        "--assembled-transcripts-support-tsv-file",
        dest="assembled_transcripts_support_tsv_file",
        type=str,
        default="",
        help="Input RNA assembly support TSV file for transcripts assembled outside exacto "
             "(columns: assembled_transcript_name, sequence, read_names)."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--stitched-transcripts-tsv-file",
        dest="stitched_transcripts_tsv_file",
        type=str,
        default="",
        required=False,
        help="Input stitched transcripts TSV file written by exacto stitch-reference-transcripts "
             "(optional). Its stitched sequences replace the transcript sequences, joined by "
             "transcript name, so the call-rna-transcript-vars TSV files must come from the "
             "realigned stitched FASTA of the same stitch run. Adds reference-stitched provenance "
             "columns to both outputs: is_reference_stitched per nucleotide, and "
             "reference_stitched_amino_acid_intervals plus num_reference_stitched_amino_acids per "
             "proteoform."
    )
    parser_optional.add_argument(
        "--start-codons",
        dest="start_codons",
        type=str,
        nargs='+',
        default=TRANSLATE_START_CODONS,
        required=False,
        help="One or more start codons. Pass multiple as: --start-codons AUG GUG CUG (default: %s)."
             % ' '.join(TRANSLATE_START_CODONS)
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=TRANSLATE_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % TRANSLATE_NUM_THREADS
    )
    parser.set_defaults(which='translate-transcripts')
    return sub_parsers


def run_cli_translate_transcripts_from_parsed_args(args) -> None:
    """
    Run Exacto 'translate-transcripts' command using parameters from parsed arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    translate_transcripts(
        assembled_transcript_support_tsv_file=args.assembled_transcripts_support_tsv_file,
        rna_consensus_tsv_file=args.rna_consensus_tsv_file,
        stitched_transcripts_tsv_file=args.stitched_transcripts_tsv_file,
        assembled_transcript_model_alignments_tsv_file=args.assembled_transcript_model_alignments_tsv_file,
        assembled_transcript_variants_tsv_file=args.assembled_transcript_variants_tsv_file,
        dna_variants_tsv_file=args.dna_variants_tsv_file,
        integrated_variants_tsv_file=args.integrated_variants_tsv_file,
        strategy=TranslationStrategy(args.strategy),
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        start_codons=args.start_codons,
        num_threads=args.num_threads
    )
