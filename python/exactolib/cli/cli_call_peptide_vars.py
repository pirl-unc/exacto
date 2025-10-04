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
and run Exacto 'call-peptide-vars' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_call_peptide_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'call-peptide-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('call-peptide-vars', help='Call peptide variants.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--fasta-file",
        dest="fasta_file",
        type=str,
        required=True,
        help="Peptides FASTA file."
    )
    parser_required.add_argument(
        "--reference-fasta-file",
        dest="reference_fasta_file",
        type=str,
        required=True,
        help="Reference peptides (e.g. GENCODE) FASTA file."
    )
    parser_required.add_argument(
        "--translations-tsv-file",
        dest="translations_tsv_file",
        type=str,
        required=True,
        help="Translations TSV file."
    )
    parser_required.add_argument(
        "--rna-variants-tsv-file",
        dest="rna_variants_tsv_file",
        type=str,
        required=True,
        help="RNA variants TSV file."
    )
    parser_required.add_argument(
        "--dna-variants-tsv-file",
        dest="dna_variants_tsv_file",
        type=str,
        required=True,
        help="DNA variants TSV file."
    )
    parser_required.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser_required.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=True,
        help="Output FASTA file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_required.add_argument(
        "--rna-bam-file",
        dest="rna_bam_file",
        type=str,
        required=False,
        default="",
        help="RNA BAM file."
    )
    parser_required.add_argument(
        "--rna-bam-bai-file",
        dest="rna_bam_bai_file",
        type=str,
        required=False,
        default="",
        help="RNA BAM.BAI file."
    )
    parser_required.add_argument(
        "--exclude-bed-file",
        dest="exclude_bed_file",
        type=str,
        default="",
        required=False,
        help="A list of regions to exclude (e.g. immunoglobulin gene regions)."
    )
    parser_optional.add_argument(
        "--k",
        dest="k",
        type=int,
        default=CALL_PEPTIDE_VARS_K,
        required=False,
        help="K-mer k size (default: %i)."
             % CALL_PEPTIDE_VARS_K
    )
    parser_optional.add_argument(
        "--min-reads",
        dest="min_reads",
        type=int,
        default=CALL_PEPTIDE_VARS_MIN_READS,
        required=False,
        help="Minimum number of translated peptide sequence (RNA) reads (default: %i)."
             % CALL_PEPTIDE_VARS_MIN_READS
    )
    parser_optional.add_argument(
        "--dna-variant-padding",
        dest="dna_variant_padding",
        type=int,
        default=CALL_PEPTIDE_VARS_DNA_VARIANT_PADDING,
        required=False,
        help="DNA variant padding (default: %i)."
             % CALL_PEPTIDE_VARS_DNA_VARIANT_PADDING
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_PEPTIDE_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_PEPTIDE_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        default=CALL_PEPTIDE_VARS_GZIP,
        required=False,
        help="If 'yes', gzip the output TSV file (default: %s)."
             % CALL_PEPTIDE_VARS_GZIP
    )

    parser.set_defaults(which='call-peptide-vars')
    return sub_parsers


def run_cli_call_peptide_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'call-peptide-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    fasta_file
                    reference_fasta_file
                    translations_tsv_file
                    output_tsv_file
                    output_fasta_file
                    k
                    min_reads
                    num_threads
                    gzip
    """
    pass
    # if args.gzip:
    #     if args.output_tsv_file.endswith('.gz'):
    #         output_tsv_file = args.output_tsv_file
    #     else:
    #         output_tsv_file = args.output_tsv_file + '.gz'
    #     if args.output_fasta_file.endswith('.gz'):
    #         output_fasta_file = args.output_fasta_file
    #     else:
    #         output_fasta_file = args.output_fasta_file + '.gz'
    # else:
    #     output_tsv_file = args.output_tsv_file
    #     output_fasta_file = args.output_fasta_file
    #
    # if args.exclude_bed_file != "":
    #     assert args.rna_bam_file != "", "If --excl-bed-file is specified, --rna-bam-file must be specified."
    #     assert args.rna_bam_bai_file != "", "If --excl-bed-file is specified, --rna-bam-bai-file must be specified."
    #
    # identify_peptide_variants(
    #     fasta_file=args.fasta_file,
    #     rna_bam_file=args.rna_bam_file,
    #     rna_bam_bai_file=args.rna_bam_bai_file,
    #     reference_fasta_file=args.reference_fasta_file,
    #     translations_tsv_file=args.translations_tsv_file,
    #     rna_variants_tsv_file=args.rna_variants_tsv_file,
    #     dna_variants_tsv_file=args.dna_variants_tsv_file,
    #     exclude_bed_file=args.exclude_bed_file,
    #     min_reads=args.min_reads,
    #     k=args.k,
    #     num_threads=args.num_threads,
    #     dna_variant_padding=args.dna_variant_padding,
    #     output_tsv_file=output_tsv_file,
    #     output_fasta_file=output_fasta_file,
    #     gzip=args.gzip
    # )
