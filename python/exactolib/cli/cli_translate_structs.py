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
and run Exacto 'translate-structs' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_translate_structs_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'translate-structs' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('translate-structs', help='Translate transcript structures to primary structures.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--rna-assembly-support-tsv-file",
        dest="rna_assembly_support_tsv_file",
        type=str,
        required=True,
        help="Input RNA assembly TSV file (columns: transcript_id, sequence, read_ids)."
    )
    parser_required.add_argument(
        "--transcript-model-structures-tsv-file",
        dest="transcript_model_structures_tsv_file",
        type=str,
        required=True,
        help="Input transcript structures TSV file."
    )
    parser_required.add_argument(
        "--rna-variants-tsv-file",
        dest="rna_variants_tsv_file",
        type=str,
        required=True,
        help="Input RNA variants TSV file."
    )
    parser_required.add_argument(
        "--dna-variants-tsv-file",
        dest="dna_variants_tsv_file",
        type=str,
        required=True,
        help="Input DNA variants TSV file."
    )
    parser_required.add_argument(
        "--integrated-variants-tsv-file",
        dest="integrated_variants_tsv_file",
        type=str,
        required=True,
        help="Input integrated variants TSV file."
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

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=TRANSLATE_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % TRANSLATE_NUM_THREADS
    )
    parser.set_defaults(which='translate-structs')
    return sub_parsers


def run_cli_translate_structs_from_parsed_args(args) -> None:
    """
    Run Exacto 'translate-structs' command using parameters from parsed arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    translate_structures(
        rna_assembly_support_tsv_file=args.rna_assembly_support_tsv_file,
        transcript_model_structures_tsv_file=args.transcript_model_structures_tsv_file,
        rna_variants_tsv_file=args.rna_variants_tsv_file,
        dna_variants_tsv_file=args.dna_variants_tsv_file,
        integrated_variants_tsv_file=args.integrated_variants_tsv_file,
        strategy=TranslationStrategy(args.strategy),
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        num_threads=args.num_threads
    )
