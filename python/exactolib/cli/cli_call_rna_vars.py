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
        "--output-ref-transcript-matches-tsv-file",
        dest="output_ref_transcript_matches_tsv_file",
        type=str,
        required=True,
        help="Output reference transcript matches TSV file."
    )
    parser_required.add_argument(
        "--output-exons-tsv-file",
        dest="output_exons_tsv_file",
        type=str,
        required=True,
        help="Output exons TSV file."
    )
    parser_required.add_argument(
        "--output-sj-tsv-file",
        dest="output_sj_tsv_file",
        type=str,
        required=True,
        help="Output splice junctions TSV file."
    )
    parser_required.add_argument(
        "--output-variants-tsv-file",
        dest="output_variants_tsv_file",
        type=str,
        required=True,
        help="Output variant calls TSV file."
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
        "--gzip",
        dest="gzip",
        type=str2bool,
        default=CALL_RNA_VARS_GZIP,
        required=False,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % CALL_RNA_VARS_GZIP
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
                    output_ref_transcript_matches_tsv_file
                    output_exons_tsv_file
                    output_sj_tsv_file
                    output_variants_tsv_file
                    num_threads
                    gzip
                    min_mapping_quality
                    min_average_base_quality
                    temp_dir
    """
    if args.gzip:
        if args.output_ref_transcript_matches_tsv_file.endswith('.gz'):
            output_ref_transcript_matches_tsv_file = args.output_ref_transcript_matches_tsv_file
        else:
            output_ref_transcript_matches_tsv_file = args.output_ref_transcript_matches_tsv_file + '.gz'
        if args.output_exons_tsv_file.endswith('.gz'):
            output_exons_tsv_file = args.output_exons_tsv_file
        else:
            output_exons_tsv_file = args.output_exons_tsv_file + '.gz'
        if args.output_sj_tsv_file.endswith('.gz'):
            output_sj_tsv_file = args.output_sj_tsv_file
        else:
            output_sj_tsv_file = args.output_sj_tsv_file + '.gz'
        if args.output_variants_tsv_file.endswith('.gz'):
            output_variants_tsv_file = args.output_variants_tsv_file
        else:
            output_variants_tsv_file = args.output_variants_tsv_file + '.gz'
    else:
        output_ref_transcript_matches_tsv_file = args.output_ref_transcript_matches_tsv_file
        output_exons_tsv_file = args.output_exons_tsv_file
        output_sj_tsv_file = args.output_sj_tsv_file
        output_variants_tsv_file = args.output_variants_tsv_file

    identify_rna_variants(
        bam_file=args.bam_file,
        bam_bai_file=args.bam_bai_file,
        reference_genome_fasta_file=args.reference_genome_fasta_file,
        gene_annotation_file=args.gene_annotation_file,
        gene_annotation_source=args.gene_annotation_source,
        output_ref_transcript_matches_tsv_file=output_ref_transcript_matches_tsv_file,
        output_exons_tsv_file=output_exons_tsv_file,
        output_sj_tsv_file=output_sj_tsv_file,
        output_variants_tsv_file=output_variants_tsv_file,
        gzip=args.gzip,
        min_mapping_quality=args.min_mapping_quality,
        min_average_base_quality=args.min_average_base_quality,
        num_threads=args.num_threads,
        temp_dir=args.temp_dir
    )
