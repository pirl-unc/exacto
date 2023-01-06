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
and run Exacto 'identify' command.
"""


import csv
import pandas as pd
import pysam

from ..constants import *
from ..main import *


def add_exacto_identify_arg_parser(sub_parsers):
    """
    Adds 'identify' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'identify',
        help='Identifies variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--nucleic_acid_type",
        dest="nucleic_acid_type",
        type=str,
        choices=NucleicAcidTypes.ALL,
        required=True,
        help="Nucleic acid type. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in NucleicAcidTypes.ALL))
    )
    parser_required.add_argument(
        "--bam_file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file."
    )
    parser_required.add_argument(
        "--fasta_file",
        dest="fasta_file",
        type=str,
        required=True,
        help="FASTA file."
    )
    # parser_required.add_argument(
    #     "--germline_resource_vcf_file",
    #     dest="germline_resource_vcf_file",
    #     type=str,
    #     required=True,
    #     help="Germline resource VCF file."
    # )
    parser_required.add_argument(
        "--output_tsv_file",
        dest="output_tsv_file",
        type=str,
        required=True,
        help="Output TSV file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num_cores",
        dest="num_cores",
        type=int,
        default=NUM_CORES,
        required=False,
        help="Number of cores to use (default: %i)."
             % NUM_CORES
    )
    parser.set_defaults(which='identify')
    return sub_parsers


def run_exacto_identify_from_parsed_args(args):
    """
    Run Exacto 'identify' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                nucleic_acid_type
                bam_file
                fasta_file
                output_tsv_file
                num_cores
    """
    if args.nucleic_acid_type == NucleicAcidTypes.RNA:
        bam_file = pysam.AlignmentFile(args.bam_file)
        df_variants = run_exacto_identify_rna_variants(
            bam_file=bam_file,
            num_cores=args.num_cores
        )
        df_variants.to_csv(args.output_tsv_file,
                           sep='\t',
                           quoting = csv.QUOTE_NONE,
                           index=False)

    # if args.variant_type == VariantTypes.SV:
    #     df_structural_variants = pd.read_csv(args.tsv_file, sep='\t')
    #     if args.annotation_source == AnnotationSources.ENSEMBL:
    #         df_structural_variants = run_exacto_annotate_genomic_structural_variants(
    #             df_structural_variants=df_structural_variants,
    #             annotation_source=args.annotation_source,
    #             df_gencode_genes=None,
    #             df_gencode_exons=None,
    #             ensembl_release=args.ensembl_release
    #         )
    #     elif args.annotation_source == AnnotationSources.GENCODE:
    #         df_gencode_genes, \
    #         df_gencode_transcripts, \
    #         df_gencode_exons = read_gencode_gtf_file(gencode_gtf_file=args.gencode_gtf_file)
    #         df_structural_variants = run_exacto_annotate_genomic_structural_variants(
    #             df_structural_variants=df_structural_variants,
    #             annotation_source=args.annotation_source,
    #             df_gencode_genes=df_gencode_genes,
    #             df_gencode_exons=df_gencode_exons,
    #             ensembl_release=None,
    #         )
    #     else:
    #         raise Exception(
    #             "Invalid value for '--annotation_source': %s. "
    #             "Allowed '--annotation_source' values are %s "
    #             % (args.annotation_source,
    #                ', '.join(f"'{item}'" for item in AnnotationSources.ALL))
    #         )
    #
    #     df_structural_variants.to_csv(args.output_tsv_file,
    #                                   sep='\t',
    #                                   index=False)
    # elif args.variant_type == VariantTypes.SNV_INDEL:
    #     df_small_variants = pd.read_csv(args.tsv_file, sep='\t')
    #     if args.annotation_source == AnnotationSources.ENSEMBL:
    #         df_small_variants = run_exacto_annotate_genomic_small_variants(
    #             df_small_variants=df_small_variants,
    #             annotation_source=args.annotation_source,
    #             df_gencode_genes=None,
    #             df_gencode_exons=None,
    #             ensembl_release=args.ensembl_release,
    #             perl_path=None,
    #             annovar_path=None,
    #             annovar_humandb_path=None,
    #             annovar_protocol=None,
    #             annovar_operation=None,
    #             annovar_genome_assembly=None,
    #             annovar_avinput_file=None,
    #             annovar_output_file=None
    #         )
    #     elif args.annotation_source == AnnotationSources.GENCODE:
    #         df_gencode_genes, \
    #         df_gencode_transcripts, \
    #         df_gencode_exons = read_gencode_gtf_file(gencode_gtf_file=args.gencode_gtf_file)
    #         df_small_variants = run_exacto_annotate_genomic_small_variants(
    #             df_small_variants=df_small_variants,
    #             annotation_source=args.annotation_source,
    #             df_gencode_genes=df_gencode_genes,
    #             df_gencode_exons=df_gencode_exons,
    #             ensembl_release=None,
    #             perl_path=None,
    #             annovar_path=None,
    #             annovar_humandb_path=None,
    #             annovar_protocol=None,
    #             annovar_operation=None,
    #             annovar_genome_assembly=None,
    #             annovar_avinput_file=None,
    #             annovar_output_file=None
    #         )
    #     elif args.annotation_source == AnnotationSources.ANNOVAR:
    #         write_annovar_avinput_file(
    #             tsv_file=args.tsv_file,
    #             output_avinput_file=args.output_avinput_file
    #         )
    #         df_small_variants = run_exacto_annotate_genomic_small_variants(
    #             df_small_variants=df_small_variants,
    #             annotation_source=args.annotation_source,
    #             df_gencode_genes=None,
    #             df_gencode_exons=None,
    #             ensembl_release=None,
    #             perl_path=args.perl_path,
    #             annovar_path=args.annovar_path,
    #             annovar_humandb_path=args.annovar_humandb_path,
    #             annovar_protocol=args.annovar_protocol,
    #             annovar_operation=args.annovar_operation,
    #             annovar_genome_assembly=args.annovar_genome_assembly,
    #             annovar_avinput_file=args.output_avinput_file,
    #             annovar_output_file=args.output_tsv_file
    #         )
    #     else:
    #         raise Exception(
    #             "Invalid value for '--annotation_source': %s. "
    #             "Allowed '--annotation_source' values are %s "
    #             % (args.annotation_source,
    #                ', '.join(f"'{item}'" for item in AnnotationSources.ALL))
    #         )
    #
    #     df_small_variants.to_csv(args.output_tsv_file,
    #                              sep='\t',
    #                              index=False)
