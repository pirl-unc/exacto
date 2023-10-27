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
and run Exacto 'sim-rna-variants' command.
"""


# import argparse
# from ..default import *
# from ..logging import get_logger
# from ..main import *
#
#
# logger = get_logger(__name__)
#
#
# def add_cli_sim_rna_variants_arg_parser(sub_parsers) -> argparse._SubParsersAction:
#     """
#     Adds 'sim-rna-variants' parser.
#
#     Parameters
#     ----------
#     sub_parsers  :   An instance of argparse.ArgumentParser subparsers.
#
#     Returns
#     -------
#     An instance of argparse.ArgumentParser subparsers.
#     """
#     parser = sub_parsers.add_parser('sim-rna-variants', help='Simulates RNA variants.')
#     parser._action_groups.pop()
#
#     # Required arguments
#     parser_required = parser.add_argument_group('required arguments')
#     parser_required.add_argument(
#         "--reference-exons-file",
#         dest="reference_exons_file",
#         required=True,
#         help="Reference exons TSV file. The following columns are expected: "
#              "'gene_id', 'transcript_id', 'exon_id', 'exon_number', 'chromosome', 'start', 'end', 'strand'."
#     )
#     parser_required.add_argument(
#         "--max-variants-per-rna",
#         dest="max_variants_per_rna",
#         type=int,
#         required=True,
#         help="References TSV file. The following columns are expected: "
#              "'gene_id', 'transcript_id', 'exon_id', 'exon_number', 'chromosome', 'start', 'end', 'strand'."
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
#         "--fraction-to-mutate",
#         dest="fraction_to_mutate",
#         type=float,
#         required=False,
#         default=SIM_RNA_VARIANTS_FRACTION_TO_MUTATE,
#         help="Fraction of transcripts (RNAs) to mutate in the reference file (default: %f)."
#              % SIM_RNA_VARIANTS_FRACTION_TO_MUTATE
#     )
#     parser_optional.add_argument(
#         "--max-variants-per-transcript",
#         dest="max_variants_per_transcript",
#         type=int,
#         required=False,
#         default=SIM_RNA_VARIANTS_MAX_VARIANTS_PER_TRANSCRIPT,
#         help="Maximum number of variants to simulate for a transcript (default: %i)." % SIM_RNA_VARIANTS_MAX_VARIANTS_PER_TRANSCRIPT
#     )
#     parser_optional.add_argument(
#         "--insertion-subtypes",
#         dest="insertion_subtypes",
#         type=str,
#         nargs='*',
#         required=False,
#         default=SIM_RNA_VARIANTS_INSERTION_SUBTYPES,
#         help="Insertion subtypes to simulate (default: %s). "
#              "Multiple values can be passed (e.g. '--insertion-subtypes INS_INTRONIC INS_DEL')."
#              % ' '.join(SIM_RNA_VARIANTS_INSERTION_SUBTYPES)
#     )
#     parser_optional.add_argument(
#         "--deletion-subtypes",
#         dest="deletion_subtypes",
#         type=str,
#         nargs='*',
#         required=False,
#         default=SIM_RNA_VARIANTS_INSERTION_SUBTYPES,
#         help="Deletion subtypes to simulate (default: %s). "
#              "Multiple values can be passed (e.g. '--deletion-subtypes DEL_EXONIC DEL_PARTIALEXONIC')."
#              % ' '.join(SIM_RNA_VARIANTS_DELETION_SUBTYPES)
#     )
#     parser.set_defaults(which='sim-rna-variants')
#     return sub_parsers
#
# def run_cli_sim_rna_variants_from_parsed_args(args):
#     """
#     Run Exacto 'sim-rna-variants' command using parameters from parsed arguments.
#
#     Parameters
#     ----------
#     args    :   An instance of argparse.ArgumentParser with the following variables:
#                 'reference_exons_file'
#                 'max_variants_per_rna'
#                 'output_tsv_file'
#                 'output_fasta_file'
#                 'fraction_to_mutate'
#                 'max_variants_per_transcript'
#                 'insertion_subtypes',
#                 'deletion_subtypes'
#     """
#     # Step 1. Load variants lists
#     logger.info("Started reading all TSV files")
#     pool = mp.Pool(processes=args.num_processes)
#     async_results = []
#     for tsv_file in args.tsv_file:
#         async_results.append(pool.apply_async(load_tsv_file_worker, args=(tsv_file,)))
#     pool.close()
#     pool.join()
#     variants_lists = [async_result.get() for async_result in async_results]
#     logger.info("Finished reading all TSV files")
#
#     # Step 2. Merge variants lists
#     logger.info("Started merging all variants into one list")
#     variants_list = run_exacto_merge_variants(
#         variants_lists=variants_lists,
#         max_neighbor_distance=args.max_neighbor_distance
#     )
#     logger.info("Finished merging all variants into one list")
#
#     # Step 3. Write to a TSV file
#     df_variants_list = variants_list.to_dataframe()
#     df_variants_list.sort_values(['variant_id'], inplace=True)
#     df_variants_list.to_csv(
#         args.output_tsv_file,
#         sep='\t',
#         index=False
#     )
