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
and run Exacto 'annotate' command.
"""


import pandas as pd
from ..constants import *
from ..main import *
from ..utilities.gencode_utils import *


def add_exacto_annotate_arg_parser(sub_parsers):
    """
    Adds 'annotate' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'annotate',
        help='Annotate variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--annotation_source",
        dest="annotation_source",
        type=str,
        choices=AnnotationSources.ALL,
        required=True,
        help="Annotation source. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in AnnotationSources.ALL))
    )
    parser_required.add_argument(
        '--variant_class',
        type=str,
        required=True,
        choices=VariantClasses.ALL,
        help="Variant class (%s). "
             "If the input TSV file is of structural variants, specify '%s'. "
             "If the input VCF file is of SNVs and INDELs, specify '%s'."
             % (', '.join(f"'{item}'" for item in VariantClasses.ALL),
                VariantClasses.SV,
                VariantClasses.SNV_INDEL)
    )
    parser_required.add_argument(
        "--tsv_file",
        dest="tsv_file",
        type=str,
        required=True,
        help="Input TSV file. "
             "If the input TSV file is of structural variants, "
             "the expected headers are: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2', 'sv_type' "
             "(DEL, INS, INV, DUP, BND or TRA). "
             "If the input TSV file is of single-nucleotide variants and INDELs, "
             "the expected headers are: 'chrom', 'pos', 'ref', 'alt'."
    )
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
        "--ensembl_release",
        dest="ensembl_release",
        type=int,
        required=False,
        help="Ensembl release version number "
             "(e.g. 75 for GRCh37 or 106 for GRCh38). "
             "This parameter must be supplied if "
             "--annotation_source is '%s'. "
             "Please make sure the specified "
             "ensembl version is installed using pyensembl."
             % AnnotationSources.ENSEMBL
    )
    parser_optional.add_argument(
        "--gencode_gtf_file",
        dest="gencode_gtf_file",
        type=str,
        required=False,
        help="GENCODE GTF file. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.GENCODE
    )
    parser_optional.add_argument(
        "--perl_path",
        dest="perl_path",
        type=str,
        required=False,
        help="Perl path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar_path",
        dest="annovar_path",
        type=str,
        required=False,
        help="ANNOVAR directory path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar_humandb_path",
        dest="annovar_humandb_path",
        type=str,
        required=False,
        help="ANNOVAR humandb/ directory path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--output_avinput_file",
        dest="output_avinput_file",
        type=str,
        required=False,
        help="Output ANNOVAR .avinput file. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar_genome_assembly",
        dest="annovar_genome_assembly",
        type=str,
        required=False,
        help="ANNOVAR genome assembly. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar_protocol",
        dest="annovar_protocol",
        type=str,
        required=False,
        default=','.join(list(ANNOVAR_PROTOCOL_OPERATION.keys())),
        help="ANNOVAR protocol (e.g. 'refGene,exac03'). "
             "This parameter must be supplied if "
             "--annotation_source is '%s' (default: '%s')."
             % (AnnotationSources.ANNOVAR, ','.join(list(ANNOVAR_PROTOCOL_OPERATION.keys())))
    )
    parser_optional.add_argument(
        "--annovar_operation",
        dest="annovar_operation",
        type=str,
        required=False,
        default=','.join(list(ANNOVAR_PROTOCOL_OPERATION.values())),
        help="ANNOVAR protocol (e.g. 'g,f'). "
             "This parameter must be supplied if "
             "--annotation_source is '%s' (default: '%s')."
             % (AnnotationSources.ANNOVAR, ','.join(list(ANNOVAR_PROTOCOL_OPERATION.values())))
    )
    parser.set_defaults(which='annotate')
    return sub_parsers


def run_exacto_annotate_from_parsed_args(args):
    """
    Run Exacto 'annotate' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                annotation_source
                variant_class
                tsv_file
                output_tsv_file
                ensembl_release
                gencode_gtf_file
    """
    if args.variant_class == VariantClasses.SV:
        df_structural_variants = pd.read_csv(args.tsv_file, sep='\t')
        if args.annotation_source == AnnotationSources.ENSEMBL:
            df_structural_variants = run_exacto_annotate_genomic_structural_variants(
                df_structural_variants=df_structural_variants,
                annotation_source=args.annotation_source,
                df_gencode_genes=None,
                df_gencode_exons=None,
                ensembl_release=args.ensembl_release
            )
        elif args.annotation_source == AnnotationSources.GENCODE:
            df_gencode_genes, \
            df_gencode_transcripts, \
            df_gencode_exons = read_gencode_gtf_file(gencode_gtf_file=args.gencode_gtf_file)
            df_structural_variants = run_exacto_annotate_genomic_structural_variants(
                df_structural_variants=df_structural_variants,
                annotation_source=args.annotation_source,
                df_gencode_genes=df_gencode_genes,
                df_gencode_exons=df_gencode_exons,
                ensembl_release=None,
            )
        else:
            raise Exception(
                "Invalid value for '--annotation_source': %s. "
                "Allowed '--annotation_source' values are %s "
                % (args.annotation_source,
                   ', '.join(f"'{item}'" for item in AnnotationSources.ALL))
            )

        df_structural_variants.to_csv(args.output_tsv_file,
                                      sep='\t',
                                      index=False)
    elif args.variant_class == VariantClasses.SNV_INDEL:
        df_small_variants = pd.read_csv(args.tsv_file, sep='\t')
        if args.annotation_source == AnnotationSources.ENSEMBL:
            df_small_variants = run_exacto_annotate_genomic_small_variants(
                df_small_variants=df_small_variants,
                annotation_source=args.annotation_source,
                df_gencode_genes=None,
                df_gencode_exons=None,
                ensembl_release=args.ensembl_release,
                perl_path=None,
                annovar_path=None,
                annovar_humandb_path=None,
                annovar_protocol=None,
                annovar_operation=None,
                annovar_genome_assembly=None,
                annovar_avinput_file=None,
                annovar_output_file=None
            )
        elif args.annotation_source == AnnotationSources.GENCODE:
            df_gencode_genes, \
            df_gencode_transcripts, \
            df_gencode_exons = read_gencode_gtf_file(gencode_gtf_file=args.gencode_gtf_file)
            df_small_variants = run_exacto_annotate_genomic_small_variants(
                df_small_variants=df_small_variants,
                annotation_source=args.annotation_source,
                df_gencode_genes=df_gencode_genes,
                df_gencode_exons=df_gencode_exons,
                ensembl_release=None,
                perl_path=None,
                annovar_path=None,
                annovar_humandb_path=None,
                annovar_protocol=None,
                annovar_operation=None,
                annovar_genome_assembly=None,
                annovar_avinput_file=None,
                annovar_output_file=None
            )
        elif args.annotation_source == AnnotationSources.ANNOVAR:
            write_annovar_avinput_file(
                tsv_file=args.tsv_file,
                output_avinput_file=args.output_avinput_file
            )
            df_small_variants = run_exacto_annotate_genomic_small_variants(
                df_small_variants=df_small_variants,
                annotation_source=args.annotation_source,
                df_gencode_genes=None,
                df_gencode_exons=None,
                ensembl_release=None,
                perl_path=args.perl_path,
                annovar_path=args.annovar_path,
                annovar_humandb_path=args.annovar_humandb_path,
                annovar_protocol=args.annovar_protocol,
                annovar_operation=args.annovar_operation,
                annovar_genome_assembly=args.annovar_genome_assembly,
                annovar_avinput_file=args.output_avinput_file,
                annovar_output_file=args.output_tsv_file
            )
        else:
            raise Exception(
                "Invalid value for '--annotation_source': %s. "
                "Allowed '--annotation_source' values are %s "
                % (args.annotation_source,
                   ', '.join(f"'{item}'" for item in AnnotationSources.ALL))
            )

        df_small_variants.to_csv(args.output_tsv_file,
                                 sep='\t',
                                 index=False)
