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
and run Exacto 'translate-seqs' command.
"""


import argparse
from ..constants import *
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_translate_seqs_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'translate-seqs' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('translate-seqs', help='Translate transcript sequences in a long-read RNA-seq FASTA or FASTQ file to peptide sequences.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')

    # Mutually exclusive group for input
    input_group = parser_required.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--fastq-file",
        dest="fastq_file",
        type=str,
        help="Input FASTQ or FASTQ.GZ file."
    )
    input_group.add_argument(
        "--fasta-file",
        dest="fasta_file",
        type=str,
        help="Input FASTA or FASTA.GZ file."
    )
    input_group.add_argument(
        "--sequence",
        dest="sequence",
        type=str,
        help="Directly input a transcript sequence."
    )

    parser_required.add_argument(
        "--strategy",
        dest="strategy",
        type=str,
        choices=[str(TranslationStrategy.LONGEST_ORF), str(TranslationStrategy.ALL_ORFS)],
        required=True,
        help="Translation strategy (default: %s). Available options: '%s'." %
             (str(TranslationStrategy.LONGEST_ORF),
              ','.join([str(TranslationStrategy.LONGEST_ORF)]))
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--output-tsv-file",
        dest="output_tsv_file",
        type=str,
        required=False,
        help="Output TSV file."
    )
    parser_optional.add_argument(
        "--output-fasta-file",
        dest="output_fasta_file",
        type=str,
        required=False,
        help="Output FASTA file."
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=TRANSLATE_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % TRANSLATE_NUM_THREADS
    )
    parser_optional.add_argument(
        "--temp-dir",
        dest="temp_dir",
        type=str,
        default="",
        required=False,
        help="Temp directory (default: TMPDIR)."
    )
    parser_optional.add_argument(
        "--gzip",
        dest="gzip",
        type=str2bool,
        default=TRANSLATE_GZIP,
        required=False,
        help="If 'yes', gzip the output TSV and FASTA file (default: %s)."
             % TRANSLATE_GZIP
    )
    parser.set_defaults(which='translate-seqs')
    return sub_parsers


def run_cli_translate_seqs_from_parsed_args(args) -> None:
    """
    Run Exacto 'translate-seqs' command using parameters from parsed arguments.

    Parameters:
        args    :   An instance of argparse.ArgumentParser with the following variables:
                    fastq_file
                    strategy
                    output_tsv_file
                    output_fasta_file
                    num_threads
                    gzip
    """
    if args.fasta_file or args.fastq_file:
        if args.output_tsv_file is None:
            raise Exception('--output-tsv-file must be specified.')
        if args.output_fasta_file is None:
            raise Exception('--output-fasta-file must be specified.')

        # Step 1. Translate
        if args.fasta_file:
            logger.info("%i reads in total in the FASTQ file." % count_reads_in_fastx_file(fastx_file=args.fasta_file))
            df_translations = translate_fasta_file(
                fasta_file=args.fasta_file,
                strategy=TranslationStrategy(args.strategy),
                num_threads=args.num_threads,
                temp_dir=args.temp_dir
            )
        elif args.fastq_file:
            logger.info("%i reads in total in the FASTQ file." % count_reads_in_fastx_file(fastx_file=args.fastq_file))
            df_translations = translate_fastq_file(
                fastq_file=args.fastq_file,
                strategy=TranslationStrategy(args.strategy),
                num_threads=args.num_threads,
                temp_dir=args.temp_dir
            )
        else:
            raise Exception("Unexpected error - either a FASTA or FASTQ file should have been specified.")
        logger.info("%i translated reads." % len(df_translations))
        logger.info('%i unique rna IDs in the translated peptides.' % len(df_translations['rna_id'].unique()))

        # Step 2. Prepare file paths with appropriate extensions
        if args.gzip:
            if args.output_tsv_file.endswith('.gz'):
                output_tsv_file = args.output_tsv_file
            else:
                output_tsv_file = args.output_tsv_file + '.gz'
            if args.output_fasta_file.endswith('.gz'):
                output_fasta_file = args.output_fasta_file
            else:
                output_fasta_file = args.output_fasta_file + '.gz'
        else:
            output_tsv_file = args.output_tsv_file
            output_fasta_file = args.output_fasta_file

        # Step 3. Output TSV file
        df_translations.to_csv(
            output_tsv_file,
            sep='\t',
            index=False,
            compression='gzip' if args.gzip else None
        )

        # Step 4. Output FASTA file
        df_peptides_unique = df_translations.loc[:, ['peptide_id', 'peptide_sequence']].drop_duplicates()
        if args.gzip:
            with pysam.BGZFile(output_fasta_file, "wb") as fasta:
                for index,row in df_peptides_unique.iterrows():
                    curr_peptide_id = str(row['peptide_id']).encode()
                    curr_peptide_sequence = str(row['peptide_sequence']).encode()
                    fasta.write(b">%s\n" % curr_peptide_id)
                    fasta.write(b"%s\n" % curr_peptide_sequence)
        else:
            with open(output_fasta_file, "w") as file:
                for index,row in df_peptides_unique.iterrows():
                    curr_peptide_id = str(row['peptide_id'])
                    curr_peptide_sequence = str(row['peptide_sequence'])
                    file.write(">%s\n" % curr_peptide_id)
                    file.write("%s\n" % curr_peptide_sequence)
        pysam.faidx(output_fasta_file, rebuild=True)
    elif args.sequence:
        translations = translate_sequence(rna_sequence=args.sequence, strategy=args.strategy)
        print('Translated peptide sequence(s):')
        print('[orf_start:orf_end] [sequence]')
        for (sequence, orf_start, orf_end) in translations:
            print('%i:%i %s' % (orf_start, orf_end, sequence))
    else:
        raise Exception("Unexpected error: one of the following should have been specified: --fasta-file, --fastq-file, or --sequence.")
