# # Licensed under the Apache License, Version 2.0 (the "License");
# # you may not use this file except in compliance with the License.
# # You may obtain a copy of the License at
# #
# #     http://www.apache.org/licenses/LICENSE-2.0
# #
# # Unless required by applicable law or agreed to in writing, software
# # distributed under the License is distributed on an "AS IS" BASIS,
# # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# # See the License for the specific language governing permissions and
# # limitations under the License.
#
#
# """
# The purpose of this python3 script is to create parser
# and run Exacto 'report' command.
# """
#
#
# import argparse
# import gzip
# import pandas as pd
# from collections import defaultdict
# from ..constants import *
# from ..default import *
# from ..logging import get_logger
# from ..main import *
# from ..utilities import *
#
#
# logger = get_logger(__name__)
#
#
# def add_cli_report_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
#     """
#     Add 'report' parser.
#
#     Parameters:
#         sub_parsers     :  argparse.ArgumentParser subparsers.
#
#     Returns:
#         sub_parsers     :   argparse.ArgumentParser subparsers
#     """
#     parser = sub_parsers.add_parser('report', help='Generate an HTML report of the DNA, RNA, and peptide variants.')
#     parser._action_groups.pop()
#
#     # Required arguments
#     parser_required = parser.add_argument_group('required arguments')
#     parser_required.add_argument(
#         "--fastq-file",
#         dest="fastq_file",
#         type=str,
#         required=True,
#         help="Input FASTQ.GZ file."
#     )
#     parser_required.add_argument(
#         "--strategy",
#         dest="strategy",
#         type=str,
#         choices=[TranslationStrategies.LONGEST_ORF],
#         required=True,
#         help="Translation strategy (default: %s). Available options: '%s'." %
#              (TranslationStrategies.LONGEST_ORF,
#               ','.join(TranslationStrategies.ALL))
#     )
#     parser_required.add_argument(
#         "--output-tsv-file",
#         dest="output_tsv_file",
#         type=str,
#         required=True,
#         help="Output TSV file."
#     )
#     parser_required.add_argument(
#         "--output-fasta-file",
#         dest="output_fasta_file",
#         type=str,
#         required=True,
#         help="Output FASTA file."
#     )
#
#     # Optional arguments
#     parser_optional = parser.add_argument_group('optional arguments')
#     parser_optional.add_argument(
#         "--num-threads",
#         dest="num_threads",
#         type=int,
#         default=TRANSLATE_NUM_THREADS,
#         required=False,
#         help="Number of threads (default: %i)."
#              % TRANSLATE_NUM_THREADS
#     )
#     parser_optional.add_argument(
#         "--temp-dir",
#         dest="temp_dir",
#         type=str,
#         default="",
#         required=False,
#         help="Temp directory (default: TMPDIR)."
#     )
#     parser_optional.add_argument(
#         "--gzip",
#         dest="gzip",
#         type=str2bool,
#         default=TRANSLATE_GZIP,
#         required=False,
#         help="If 'yes', gzip the output TSV and FASTA file (default: %s)."
#              % TRANSLATE_GZIP
#     )
#     parser.set_defaults(which='translate')
#     return sub_parsers
#
#
# def run_cli_translate_vars_from_parsed_args(args) -> None:
#     """
#     Run Exacto 'translate' command using parameters from parsed arguments.
#
#     Parameters:
#         args    :   An instance of argparse.ArgumentParser with the following variables:
#                     fastq_file
#                     strategy
#                     output_tsv_file
#                     output_fasta_file
#                     num_threads
#                     gzip
#     """
#     # Step 1. Count the number of reads in the FASTQ file
#     num_reads = count_reads_in_fastq(fastq_file=args.fastq_file)
#     logger.info("%i reads in total in the FASTQ file." % num_reads)
#
#     # Step 2. Translate
#     df_translations = translate(
#         fastq_file=args.fastq_file,
#         strategy=args.strategy,
#         num_threads=args.num_threads,
#         temp_dir=args.temp_dir
#     )
#     logger.info("%i translated reads." % len(df_translations))
#     logger.info('%i unique read IDs in the translated peptides.' % len(df_translations['read_name'].unique()))
#
#     # Step 3. Prepare file paths with appropriate extensions
#     if args.gzip:
#         if args.output_tsv_file.endswith('.gz'):
#             output_tsv_file = args.output_tsv_file
#         else:
#             output_tsv_file = args.output_tsv_file + '.gz'
#         if args.output_fasta_file.endswith('.gz'):
#             output_fasta_file = args.output_fasta_file
#         else:
#             output_fasta_file = args.output_fasta_file + '.gz'
#     else:
#         output_tsv_file = args.output_tsv_file
#         output_fasta_file = args.output_fasta_file
#
#     # Step 4. Output TSV file
#     df_translations.to_csv(
#         output_tsv_file,
#         sep='\t',
#         index=False,
#         compression='gzip' if args.gzip else None
#     )
#
#     # Step 5. Output FASTA file
#     df_peptides_unique = df_translations.loc[:, ['peptide_id', 'peptide_sequence']].drop_duplicates()
#     if args.gzip:
#         with pysam.BGZFile(output_fasta_file, "wb") as fasta:
#             for index,row in df_peptides_unique.iterrows():
#                 curr_peptide_id = str(row['peptide_id']).encode()
#                 curr_peptide_sequence = str(row['peptide_sequence']).encode()
#                 fasta.write(b">%s\n" % curr_peptide_id)
#                 fasta.write(b"%s\n" % curr_peptide_sequence)
#     else:
#         with open(output_fasta_file, "w") as file:
#             for index,row in df_peptides_unique.iterrows():
#                 curr_peptide_id = str(row['peptide_id'])
#                 curr_peptide_sequence = str(row['peptide_sequence'])
#                 file.write(">%s\n" % curr_peptide_id)
#                 file.write("%s\n" % curr_peptide_sequence)
#     pysam.faidx(output_fasta_file, rebuild=True)
