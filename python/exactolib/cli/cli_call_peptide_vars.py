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
    parser = sub_parsers.add_parser(
        'call-peptide-vars',
        help='Call peptide variants from primary structures.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--primary-structures-tsv-file",
        dest="primary_structures_tsv_file",
        type=str,
        required=True,
        help="Input primary structures TSV file."
    )
    parser_required.add_argument(
        "--reference-fasta-file",
        dest="reference_fasta_file",
        type=str,
        required=True,
        help="Reference proteome FASTA file."
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
    parser_optional.add_argument(
        "--min-k",
        dest="min_k",
        type=int,
        default=CALL_PEPTIDE_VARS_MIN_K,
        required=False,
        help="Minimum peptide length (default: %i)."
             % CALL_PEPTIDE_VARS_MIN_K
    )
    parser_optional.add_argument(
        "--max-k",
        dest="max_k",
        type=int,
        default=CALL_PEPTIDE_VARS_MAX_K,
        required=False,
        help="Maximum peptide length (default: %i)."
             % CALL_PEPTIDE_VARS_MAX_K
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

    parser.set_defaults(which='call-peptide-vars')
    return sub_parsers


def run_cli_call_peptide_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'call-peptide-vars' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    primary_structures_tsv_file
                    reference_fasta_file
                    output_tsv_file
                    min_k
                    max_k
                    num_threads
    """
    logger.info("Primary structure TSV: %s" % args.primary_structures_tsv_file)
    logger.info("Reference proteome FASTA: %s" % args.reference_fasta_file)
    logger.info("Peptide length range: [%i, %i]" % (args.min_k, args.max_k))
    logger.info("Number of threads: %i" % args.num_threads)

    df = identify_peptide_variants(
        primary_structures_tsv_file=args.primary_structures_tsv_file,
        reference_proteome_fasta_file=args.reference_fasta_file,
        min_k=args.min_k,
        max_k=args.max_k,
        num_processes=args.num_threads
    )

    df.to_csv(args.output_tsv_file, sep='\t', index=False)

    # One FASTA record per unique mutant peptide. Header lists every primary
    # structure the k-mer was observed in.
    with open(args.output_fasta_file, 'w') as f:
        for mutant_peptide_id, group in df.groupby('mutant_peptide_id'):
            sequence = group.iloc[0]['mutant_peptide_sequence']
            primary_structure_ids = ','.join(
                str(psid) for psid in sorted(group['primary_structure_id'].unique())
            )
            f.write('>mutant_peptide_id=%s primary_structure_ids=%s\n%s\n' % (
                mutant_peptide_id, primary_structure_ids, sequence
            ))

    logger.info("Wrote %i mutant peptide(s) to %s" % (len(df), args.output_tsv_file))
    logger.info("Wrote %i unique sequence(s) to %s" % (df['mutant_peptide_id'].nunique(), args.output_fasta_file))
