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
and run Exacto 'call-rna-vars' command.
"""


import argparse
import os

from ..constants import GeneAnnotationSources
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_call_rna_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'call-rna-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('call-rna-vars', help='Call RNA variants in a long-read RNA-seq BAM file.')
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
        "--reference-genome-fasta-file",
        dest="reference_genome_fasta_file",
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    parser_required.add_argument(
        "--gene-annotation-file",
        dest="gene_annotation_file",
        type=str,
        required=True,
        help="Reference gene annotation file."
    )
    parser_required.add_argument(
        "--gene-annotation-source",
        dest="gene_annotation_source",
        type=str,
        required=True,
        help="Reference gene annotation source (choices: %s)." %
             ','.join(GeneAnnotationSources.ALL)
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
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_RNA_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_RNA_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % CALL_RNA_VARS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--reference-transcript-scoring-method",
        dest="reference_transcript_scoring_method",
        type=str,
        default=CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SCORING_METHOD,
        required=False,
        help="Reference transcript scoring method (default: %s)."
             % CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SCORING_METHOD
    )
    parser_optional.add_argument(
        "--min-average-base-quality",
        dest="min_average_base_quality",
        type=float,
        default=CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY,
        required=False,
        help="Minimum average base quality (default: %f)."
             % CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY
    )
    parser_optional.add_argument(
        "--temp-dir",
        dest="temp_dir",
        type=str,
        default="",
        required=False,
        help="Temp directory (default: TMPDIR)."
    )
    parser.set_defaults(which='call-rna-vars')
    return sub_parsers


def run_cli_call_rna_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'call-rna-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    bam_file
                    bam_bai_file
                    reference_genome_fasta_file
                    gene_annotation_file
                    gene_annotation_source
                    output_dir
                    output_prefix
                    num_threads
                    min_mapping_quality
                    min_average_base_quality
                    temp_dir
    """
    os.makedirs(args.output_dir, exist_ok=True)
    identify_rna_variants(
        bam_file=args.bam_file,
        bam_bai_file=args.bam_bai_file,
        reference_genome_fasta_file=args.reference_genome_fasta_file,
        gene_annotation_file=args.gene_annotation_file,
        gene_annotation_source=args.gene_annotation_source,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        reference_transcript_scoring_method=args.reference_transcript_scoring_method,
        min_mapping_quality=args.min_mapping_quality,
        min_average_base_quality=args.min_average_base_quality,
        num_threads=args.num_threads,
        temp_dir=args.temp_dir
    )
