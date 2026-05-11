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
and run Exacto 'integrate-vars' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_integrate_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'integrate-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('integrate-vars', help='Integrate DNA and RNA variants.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--annotated-dna-vars-tsv-file",
        dest="annotated_dna_vars_tsv_file",
        type=str,
        required=True,
        help="Annotated DNA variant callset TSV file."
    )
    parser_required.add_argument(
        "--rna-vars-tsv-file",
        dest="rna_vars_tsv_file",
        type=str,
        required=True,
        help="RNA variant callset TSV file.."
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
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=INTEGRATE_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % INTEGRATE_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--max-exon-offset",
        dest="max_exon_offset",
        type=int,
        default=INTEGRATE_VARS_MAX_EXON_OFFSET,
        required=False,
        help="Maximum exon offset (default: %i)."
             % INTEGRATE_VARS_MAX_EXON_OFFSET
    )
    parser_optional.add_argument(
        "--max-transcript-boundary-offset",
        dest="max_transcript_boundary_offset",
        type=int,
        default=INTEGRATE_VARS_MAX_TRANSCRIPT_BOUNDARY_OFFSET,
        required=False,
        help="Maximum transcript boundary offset (default: %i)."
             % INTEGRATE_VARS_MAX_TRANSCRIPT_BOUNDARY_OFFSET
    )
    parser_optional.add_argument(
        "--max-intergenic-distance",
        dest="max_intergenic_distance",
        type=int,
        default=INTEGRATE_VARS_MAX_INTERGENIC_DISTANCE,
        required=False,
        help="Maximum intergenic distance (default: %i)."
             % INTEGRATE_VARS_MAX_INTERGENIC_DISTANCE
    )
    parser.set_defaults(which='integrate-vars')
    return sub_parsers


def run_cli_integrate_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'integrate-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    annotated_dna_vars_tsv_file
                    rna_vars_tsv_file
                    reference_gene_annotation_file
                    reference_gene_annotation_source
                    output_tsv_file
                    num_threads
                    max_exon_offset
                    max_transcript_boundary_offset
                    max_intergenic_distance
                    temp_dir
    """
    integrate_variants(
        dna_variant_call_annotation_set_tsv_file=args.annotated_dna_vars_tsv_file,
        rna_variant_call_set_tsv_file=args.rna_vars_tsv_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_source=args.reference_gene_annotation_source,
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        output_tsv_file=args.output_tsv_file,
        max_exon_offset=args.max_exon_offset,
        max_transcript_boundary_offset=args.max_transcript_boundary_offset,
        max_intergenic_distance=args.max_intergenic_distance,
        num_threads=args.num_threads
    )

