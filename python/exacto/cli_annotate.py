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


import argparse
from exacto.main import *
from exacto.gencode import *


def add_cli_annotate_arg_parser(sub_parsers) -> argparse._SubParsersAction:
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
        "--tsv-file",
        dest="tsv_file",
        type=str,
        required=True,
        help="Input TSV file. "
             "The expected headers are: "
             "'chr_1', 'pos_1', 'chr_2', 'pos_2'."
    )
    parser_required.add_argument(
        "--annotation-source",
        dest="annotation_source",
        type=str,
        choices=AnnotationSources.ALL,
        required=True,
        help="Annotation source. Allowed options: %s."
             % (', '.join(f"'{item}'" for item in AnnotationSources.ALL))
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
        "--ensembl-release",
        dest="ensembl_release",
        type=int,
        required=False,
        help="Ensembl release version number "
             "(e.g. 75 for GRCh37 or 106 for GRCh38). "
             "This parameter must be supplied if "
             "--annotation-source is '%s'. "
             "Please make sure the specified "
             "ensembl version is installed using pyensembl."
             % AnnotationSources.ENSEMBL
    )
    parser_optional.add_argument(
        "--ensembl-species",
        dest="ensembl_species",
        type=str,
        required=False,
        help="Ensembl species "
             "(e.g. 'human' or 'mouse'). "
             "This parameter must be supplied if "
             "--annotation-source is '%s'. "
             "Please make sure the specified "
             "ensembl species is installed using pyensembl."
             % AnnotationSources.ENSEMBL
    )
    parser_optional.add_argument(
        "--gencode-gtf-file",
        dest="gencode_gtf_file",
        type=str,
        required=False,
        help="GENCODE GTF file. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.GENCODE
    )
    parser_optional.add_argument(
        "--perl-path",
        dest="perl_path",
        type=str,
        required=False,
        help="Perl path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar-path",
        dest="annovar_path",
        type=str,
        required=False,
        help="ANNOVAR directory path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar-humandb-path",
        dest="annovar_humandb_path",
        type=str,
        required=False,
        help="ANNOVAR humandb/ directory path. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--output-avinput-file",
        dest="output_avinput_file",
        type=str,
        required=False,
        help="Output ANNOVAR .avinput file. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar-genome-assembly",
        dest="annovar_genome_assembly",
        type=str,
        required=False,
        help="ANNOVAR genome assembly. "
             "This parameter must be supplied if "
             "--annotation_source is '%s'."
             % AnnotationSources.ANNOVAR
    )
    parser_optional.add_argument(
        "--annovar-protocol",
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
        "--annovar-operation",
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


def run_cli_annotate_from_parsed_args(args) -> None:
    """
    Run Exacto 'annotate' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser
                with the following variables:
                tsv_file
                annotation_source
                output_tsv_file
                ensembl_release
                ensembl_species
                gencode_gtf_file
                perl_path
                annovar_path
                annovar_humandb_path
                output_avinput_file
                annovar_genome_assembly
                annovar_protocol
                annovar_operation
    """
    # Step 1. Load variants
    variants_list = VariantsList.read_tsv_file(tsv_file=args.tsv_file)

    # Step 2. Load annotation data
    if args.annotation_source == AnnotationSources.ENSEMBL:
        annotation = Ensembl(
            source=AnnotationSources.ENSEMBL,
            release=args.ensembl_release,
            species=args.ensembl_species
        )
    if args.annotation_source == AnnotationSources.GENCODE:
        annotation = Gencode(source=AnnotationSources.GENCODE)
        annotation.read_comprehensive_gene_annotation_gtf_file(
            gtf_file=args.gencode_gtf_file
        )

    # Step 3. Annotate variants
    variants_list = run_exacto_annotate(
        variants_list=variants_list,
        annotation=annotation
    )

    # Step 4. Write to a TSV file
    variants_list.to_dataframe().to_csv(
        args.output_tsv_file,
        sep='\t',
        index=False
    )
