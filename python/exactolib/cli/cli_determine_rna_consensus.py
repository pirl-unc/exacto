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
and run Exacto 'determine-rna-consensus' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_determine_rna_consensus_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'determine-rna-consensus' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser(
        'determine-rna-consensus',
        help="Determine an RNA consensus sequence for each RNA read cluster.")
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--tsv-file",
        dest="tsv_file",
        type=str,
        required=True,
        help="Input TSV file. Expected columns: 'cluster_id', 'read_name'. This is one of the outputs from exacto cluster-rna-reads."
    )
    parser_required.add_argument(
        "--fastq-file",
        dest="fastq_file",
        type=str,
        required=True,
        help="Input FASTQ file."
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
        "--match-score",
        dest="match_score",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_SCORE_MATCH,
        required=False,
        help="Match score for partial order alignment, added on a match. abPOA scores in "
             "positive magnitudes; a negative value is read as its magnitude "
             "(default: %i)." % DETERMINE_RNA_CONSENSUS_SCORE_MATCH
    )
    parser_optional.add_argument(
        "--mismatch-score",
        dest="mismatch_score",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_SCORE_MISMATCH,
        required=False,
        help="Mismatch penalty for partial order alignment, a positive magnitude that is "
             "subtracted (default: %i)." % DETERMINE_RNA_CONSENSUS_SCORE_MISMATCH
    )
    parser_optional.add_argument(
        "--gap-open-score",
        dest="gap_open_score",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_SCORE_GAP_OPEN,
        required=False,
        help="Gap-open penalty for partial order alignment, a positive magnitude; a gap of "
             "length L costs gap_open + gap_extend * L. Zero selects abPOA's linear gap model; "
             "a negative value is refused (default: %i)." % DETERMINE_RNA_CONSENSUS_SCORE_GAP_OPEN
    )
    parser_optional.add_argument(
        "--gap-extend-score",
        dest="gap_extend_score",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_SCORE_GAP_EXTEND,
        required=False,
        help="Gap-extend penalty for partial order alignment, a positive magnitude charged per "
             "gapped base; a value below 1 is refused (default: %i)."
             % DETERMINE_RNA_CONSENSUS_SCORE_GAP_EXTEND
    )
    parser_optional.add_argument(
        "--max-reads-per-cluster",
        dest="max_reads_per_cluster",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_MAX_READS_PER_CLUSTER,
        required=False,
        help="Maximum number of reads aligned for one cluster. A cluster of more reads gets "
             "its consensus sequence from this many of its reads, drawn at random. 0 aligns "
             "all the reads of every cluster (default: %i)."
             % DETERMINE_RNA_CONSENSUS_MAX_READS_PER_CLUSTER
    )
    parser_optional.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_SEED,
        required=False,
        help="Seed for the reads drawn from a cluster of more than --max-reads-per-cluster "
             "reads (default: %i)." % DETERMINE_RNA_CONSENSUS_SEED
    )
    parser_optional.add_argument(
        "--orientation-kmer-size",
        dest="orientation_kmer_size",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_ORIENTATION_KMER_SIZE,
        required=False,
        help="Length of the k-mers by which the reads of a cluster are turned to the "
             "orientation of most of them before alignment, at most 32. 0 aligns the reads "
             "as they are in the FASTQ file (default: %i)."
             % DETERMINE_RNA_CONSENSUS_ORIENTATION_KMER_SIZE
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=DETERMINE_RNA_CONSENSUS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % DETERMINE_RNA_CONSENSUS_NUM_THREADS
    )

    parser.set_defaults(which='determine-rna-consensus')
    return sub_parsers


def run_cli_determine_rna_consensus_from_parsed_args(args) -> None:
    """
    Run Exacto 'determine-rna-consensus' command using parameters from parsed arguments.
    """
    determine_rna_consensus(
        tsv_file=args.tsv_file,
        fastq_file=args.fastq_file,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_open_score=args.gap_open_score,
        gap_extend_score=args.gap_extend_score,
        max_reads_per_cluster=args.max_reads_per_cluster,
        seed=args.seed,
        orientation_kmer_size=args.orientation_kmer_size,
        num_threads=args.num_threads
    )

