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
and run Exacto 'remove-unspliced-rnas' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_remove_unspliced_rnas_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'remove-unspliced-rnas' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser(
        'remove-unspliced-rnas',
        help='Remove unspliced RNAs.',
        description="Remove unspliced RNAs from a coordinate-sorted BAM file of RNA reads (or assembled transcripts). A read is kept when one of its primary or supplementary alignments reaches the minimum mapping quality and either splices (has an N CIGAR operation) or overlaps a single-exon reference transcript. Secondary alignments are ignored and are not written. The output BAM holds the kept reads' primary and supplementary alignments, in the order of the input; an output file that is one of the input files is refused. The single-exon transcripts are those that pass the four annotation filters. At the default filter (protein_coding, levels 1 and 2) every read on chrM is removed: the mitochondrial protein-coding genes are GENCODE level 3 and their mRNAs are unspliced. Exacto translates with the standard genetic code only, so mitochondrial proteins are out of scope."
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input coordinate-sorted BAM file of RNA reads (or assembled transcripts)."
    )
    parser_required.add_argument(
        "--bam-bai-file",
        dest="bam_bai_file",
        type=str,
        required=True,
        help="Input BAM.BAI file."
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
        "--output-bam-file",
        dest="output_bam_file",
        type=str,
        required=True,
        help="Output BAM file (must not be the input BAM file)."
    )
    parser_required.add_argument(
        "--output-bam-bai-file",
        dest="output_bam_bai_file",
        type=str,
        required=True,
        help="Output BAM.BAI file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=REMOVE_UNSPLICED_RNAS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % REMOVE_UNSPLICED_RNAS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=REMOVE_UNSPLICED_RNAS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % REMOVE_UNSPLICED_RNAS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--gene-types",
        dest="gene_types",
        type=str,
        nargs="*",
        default=['protein_coding'],
        required=False,
        help="Reference gene types to include in annotation; give the flag with no values for no filter (default: ['protein_coding'])."
    )
    parser_optional.add_argument(
        "--gene-levels",
        dest="gene_levels",
        type=int,
        nargs="*",
        default=[1,2],
        required=False,
        help="Reference gene levels to include in annotation; give the flag with no values for no filter (default: [1,2])."
    )
    parser_optional.add_argument(
        "--transcript-types",
        dest="transcript_types",
        type=str,
        nargs="*",
        default=['protein_coding'],
        required=False,
        help="Reference transcript types to include in annotation; give the flag with no values for no filter (default: ['protein_coding'])."
    )
    parser_optional.add_argument(
        "--transcript-levels",
        dest="transcript_levels",
        type=int,
        nargs="*",
        default=[1,2],
        required=False,
        help="Reference transcript levels to include in annotation; give the flag with no values for no filter (default: [1,2])."
    )
    parser.set_defaults(which='remove-unspliced-rnas')
    return sub_parsers


def run_cli_remove_unspliced_rnas_from_parsed_args(args) -> None:
    """
    Run Exacto 'remove-unspliced-rnas' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    bam_file
                    reference_gene_annotation_file
                    reference_gene_annotation_source
                    output_bam_file
                    output_tsv_file
                    num_threads
                    min_mapping_quality
    """
    remove_unspliced_rnas(
        bam_file=args.bam_file,
        bam_bai_file=args.bam_bai_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_source=GeneAnnotationSource(args.reference_gene_annotation_source),
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        gene_types=args.gene_types,
        gene_levels=args.gene_levels,
        transcript_types=args.transcript_types,
        transcript_levels=args.transcript_levels,
        output_bam_file=args.output_bam_file,
        output_bam_bai_file=args.output_bam_bai_file,
        num_threads=args.num_threads,
        min_mapping_quality=args.min_mapping_quality
    )
