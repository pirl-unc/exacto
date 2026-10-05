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
The purpose of this python3 script is to implement the primary Exacto command.
"""


import exactolib
from .cli_annotate_vars import *
from .cli_build_genome_var_graph import *
from .cli_build_transcriptome_var_graph import *
from .cli_call_germline_dna_vars import *
from .cli_call_somatic_dna_vars import *
from .cli_call_rna_transcript_vars import *
from .cli_call_peptide_vars import *
from .cli_cluster_rna_reads import *
from .cli_correct_rna_reads import *
from .cli_determine_rna_consensus import *
from .cli_integrate_vars import *
from .cli_quantify_rna_abundances import *
from .cli_remove_unspliced_rnas import *
from .cli_stitch_reference_transcripts import *
from .cli_translate_seqs import *
from .cli_translate_transcripts import *
from ..logging import get_logger


logger = get_logger(__name__)


def init_arg_parser():
    """
    Initialize the input argument parser.

    Returns:
        argparse.ArgumentParser object
        argparse.ArgumentParser subparsers object
    """
    arg_parser = argparse.ArgumentParser(
        description="Exacto"
    )
    arg_parser.add_argument(
        '--version', '-v',
        action='version',
        version='%(prog)s version ' + str(exactolib.__version__)
    )
    sub_parsers = arg_parser.add_subparsers(help='Exacto sub-commands.')
    return arg_parser, sub_parsers


def run():
    # Step 1. Initialize argument parser
    arg_parser, sub_parsers = init_arg_parser()
    sub_parsers = add_cli_annotate_vars_arg_parser(sub_parsers=sub_parsers)                     # annotate-vars
    sub_parsers = add_cli_build_genome_var_graph_arg_parser(sub_parsers=sub_parsers)            # build-genome-var-graph
    sub_parsers = add_cli_build_transcriptome_var_graph_arg_parser(sub_parsers=sub_parsers)     # build-transcriptome-var-graph
    sub_parsers = add_cli_call_germline_dna_vars_arg_parser(sub_parsers=sub_parsers)            # call-germline-dna-vars
    sub_parsers = add_cli_call_somatic_dna_vars_arg_parser(sub_parsers=sub_parsers)             # call-somatic-dna-vars
    sub_parsers = add_cli_call_rna_transcript_vars_arg_parser(sub_parsers=sub_parsers)          # call-rna-transcript-vars
    sub_parsers = add_cli_call_peptide_vars_arg_parser(sub_parsers=sub_parsers)                 # call-peptide-vars
    sub_parsers = add_cli_cluster_rna_reads_arg_parser(sub_parsers=sub_parsers)                 # cluster-rna-reads
    sub_parsers = add_cli_correct_rna_reads_arg_parser(sub_parsers=sub_parsers)                 # correct-rna-reads
    sub_parsers = add_cli_determine_rna_consensus_arg_parser(sub_parsers=sub_parsers)           # determine-rna-consensus
    sub_parsers = add_cli_integrate_vars_arg_parser(sub_parsers=sub_parsers)                    # integrate-vars
    sub_parsers = add_cli_quantify_rna_abundances_arg_parser(sub_parsers=sub_parsers)           # quantify-rna-abundances
    sub_parsers = add_cli_remove_unspliced_rnas_arg_parser(sub_parsers=sub_parsers)             # remove-unspliced-rnas
    sub_parsers = add_cli_stitch_reference_transcripts_arg_parser(sub_parsers=sub_parsers)      # stitch-reference-transcripts
    sub_parsers = add_cli_translate_seqs_arg_parser(sub_parsers=sub_parsers)                    # translate-seqs
    sub_parsers = add_cli_translate_transcripts_arg_parser(sub_parsers=sub_parsers)             # translate-transcripts
    args = arg_parser.parse_args()

    # Step 2. Execute function based on CLI arguments
    if args.which == 'annotate-vars':
        run_cli_annotate_vars_from_parsed_args(args=args)
    elif args.which == 'build-genome-var-graph':
        run_cli_build_genome_var_graph_from_parsed_args(args=args)
    elif args.which == 'build-transcriptome-var-graph':
        run_cli_build_transcriptome_var_graph_from_parsed_args(args=args)
    elif args.which == 'call-germline-dna-vars':
        run_cli_call_germline_dna_vars_from_parsed_args(args=args)
    elif args.which == 'call-somatic-dna-vars':
        run_cli_call_somatic_dna_vars_from_parsed_args(args=args)
    elif args.which == 'call-rna-transcript-vars':
        run_cli_call_rna_transcript_vars_from_parsed_args(args=args)
    elif args.which == 'call-peptide-vars':
        run_cli_call_peptide_vars_from_parsed_args(args=args)
    elif args.which == 'cluster-rna-reads':
        run_cli_cluster_rna_reads_vars_from_parsed_args(args=args)
    elif args.which == 'correct-rna-reads':
        run_cli_correct_rna_reads_from_parsed_args(args=args)
    elif args.which == 'determine-rna-consensus':
        run_cli_determine_rna_consensus_from_parsed_args(args=args)
    elif args.which == 'integrate-vars':
        run_cli_integrate_vars_from_parsed_args(args=args)
    elif args.which == 'quantify-rna-abundances':
        run_cli_quantify_rna_abundances_vars_from_parsed_args(args=args)
    elif args.which == 'remove-unspliced-rnas':
        run_cli_remove_unspliced_rnas_from_parsed_args(args=args)
    elif args.which == 'stitch-reference-transcripts':
        run_cli_stitch_reference_transcripts_from_parsed_args(args=args)
    elif args.which == 'translate-seqs':
        run_cli_translate_seqs_from_parsed_args(args=args)
    elif args.which == 'translate-transcripts':
        run_cli_translate_transcripts_from_parsed_args(args=args)
    else:
        raise Exception("Invalid command: %s" % args.which)
