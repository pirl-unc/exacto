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
and run Exacto 'sim-variants' command.
"""


import argparse
import pysam
from .main import *
from .gencode import *
from .default_parameters import *


logger = get_logger(__name__)


def add_exacto_simulate_variants_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'sim-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'sim-variants',
        help='Simulate DNA and RNA variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        '--reference_genome_fasta_file',
        dest='reference_genome_fasta_file',
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    parser_required.add_argument(
        '--reference_transcripts_gtf_file',
        dest='reference_transcripts_gtf_file',
        type=str,
        required=True,
        help="Reference transcripts GTF file (recommended: GENCODE comprehensive gene annotation GTF file)."
    )
    parser_required.add_argument(
        '--output_tsv_file',
        dest='output_tsv_file',
        type=str,
        required=True,
        help="Output TSV file."
    )
    parser_required.add_argument(
        '--sample_id',
        dest='sample_id',
        type=str,
        required=True,
        help="Sample ID."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--target_regions_tsv_file",
        dest="target_regions_tsv_file",
        type=str,
        required=False,
        help="If this parameter is not specified, all transcripts specified in "
             "--reference_transcripts_gtf_file will be subject to variant simulation. "
             "If specified, Transcripts within the genomic regions specified in this "
             "file will be simulated for variants. Expected headers: 'chrom', 'start', 'end'"
    )
    parser_optional.add_argument(
        '--transcript_types',
        dest='transcript_types',
        nargs='+',
        required=False,
        help="Transcript types to simulate variants (e.g. 'protein_coding')."
    )
    parser_optional.add_argument(
        "--num_snv",
        dest="num_snv",
        type=int,
        default=SIMULATE_VARIANTS_NUM_SNV,
        required=False,
        help="Number of SNVs to simulate (default: %i)."
             % SIMULATE_VARIANTS_NUM_SNV
    )
    parser_optional.add_argument(
        "--num_insertion",
        dest="num_insertion",
        type=int,
        default=SIMULATE_VARIANTS_NUM_INSERTION,
        required=False,
        help="Number of insertions to simulate (default: %i)."
             % SIMULATE_VARIANTS_NUM_INSERTION
    )
    parser_optional.add_argument(
        "--insertion_size_mean",
        dest="insertion_size_mean",
        type=int,
        default=SIMULATE_VARIANTS_INSERTION_SIZE_MEAN,
        required=False,
        help="Insertion size mean (default: %i)."
             % SIMULATE_VARIANTS_INSERTION_SIZE_MEAN
    )
    parser_optional.add_argument(
        "--insertion_size_stdev",
        dest="insertion_size_stdev",
        type=int,
        default=SIMULATE_VARIANTS_INSERTION_SIZE_STDEV,
        required=False,
        help="Insertion size standard deviation (default: %i)."
             % SIMULATE_VARIANTS_INSERTION_SIZE_STDEV
    )
    parser_optional.add_argument(
        "--num_deletion",
        dest="num_deletion",
        type=int,
        default=SIMULATE_VARIANTS_NUM_DELETION,
        required=False,
        help="Number of deletions to simulate (default: %i)."
             % SIMULATE_VARIANTS_NUM_DELETION
    )
    parser_optional.add_argument(
        "--deletion_size_mean",
        dest="deletion_size_mean",
        type=int,
        default=SIMULATE_VARIANTS_DELETION_MEAN,
        required=False,
        help="Deletion size mean (default: %i)."
             % SIMULATE_VARIANTS_DELETION_MEAN
    )
    parser_optional.add_argument(
        "--deletion_size_stdev",
        dest="deletion_size_stdev",
        type=int,
        default=SIMULATE_VARIANTS_DELETION_STDEV,
        required=False,
        help="Deletion size standard deviation (default: %i)."
             % SIMULATE_VARIANTS_DELETION_STDEV
    )
    parser_optional.add_argument(
        "--enforce_infinite_sites_model",
        dest="enforce_infinite_sites_model",
        type=bool,
        default=SIMULATE_VARIANTS_ENFORCE_INFINITE_SITES_MODEL,
        required=False,
        help="If true, the simulation enforces infinite sites model (default: %r)."
             % SIMULATE_VARIANTS_ENFORCE_INFINITE_SITES_MODEL
    )
    parser.set_defaults(which='sim-variants')
    return sub_parsers


def run_exacto_sim_rna_variants_from_parsed_args(args) -> None:
    """
    Run Exacto 'sim-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                reference_genome_fasta_file
                reference_transcripts_gtf_file
                output_tsv_file
                sample_id
                target_regions_tsv_file
                transcript_types
                num_snv
                num_insertion
                num_deletion
                enforce_infinite_sites_model
    """
    logger.info("Started running exacto 'sim-variants' command.")

    # Step 1. Load input data
    logger.info("Loading reference files (reference genome FASTA and transcripts GTF files).")
    genome_fasta = pysam.FastaFile(args.reference_genome_fasta_file)
    gencode = Gencode()
    gencode.comprehensive_gene_annotation_gtf_file(gtf_file=args.reference_transcripts_gtf_file)
    if args.target_regions_tsv_file:
        logger.info("Loading target regions TSV file.")
        df_target_regions = pd.read_csv(args.target_regions_tsv_file, sep='\t')
    else:
        df_target_regions = None

    # Step 2. Simulate variants


    #
    # df_transcripts = df_transcripts.loc[df_transcripts['transcript_type'].isin(args.transcript_types),:]
    # for curr_sample_idx in range(0, args.num_samples):
    #     curr_sample_id = args.sample_id_prefix + '-' + str(curr_sample_idx + 1).zfill(4)
    #     df_rna_variants, variant_transcript_sequences = run_exacto_simulate_rna_variants(
    #         genome_fasta=genome_fasta,
    #         df_genes=df_genes,
    #         df_transcripts=df_transcripts,
    #         df_exons=df_exons,
    #         df_target_regions=df_target_regions,
    #         df_herv_regions=df_herv_regions,
    #         num_snv=args.num_snv,
    #         num_insertion=args.num_insertion,
    #         num_deletion=args.num_deletion,
    #         num_fusion=args.num_fusion,
    #         num_inversion=args.num_inversion,
    #         num_herv=args.num_herv,
    #         insertion_size_mean=args.insertion_size_mean,
    #         insertion_size_stdev=args.insertion_size_stdev,
    #         deletion_size_mean=args.deletion_size_mean,
    #         deletion_size_stdev=args.deletion_size_stdev,
    #         herv_solo_ltr_proportion=args.herv_solo_ltr_proportion,
    #         herv_truncated_proportion=args.herv_truncated_proportion,
    #         herv_chimeric_proportion=args.herv_chimeric_proportion,
    #         herv_chimeric_max_neighboring_distance=args.herv_chimeric_max_neighboring_distance,
    #         herv_full_length_proportion=args.herv_full_length_proportion,
    #         infinite_sites_assumption=args.infinite_sites_assumption
    #     )
    #
    #     # Save to files
    #     logger.info("Started writing simulated RNA variants files [%i/%i]."
    #                 % (curr_sample_idx + 1, args.num_samples))
    #     df_rna_variants.to_csv(args.output_dir + curr_sample_id + '_rna_variants.tsv',
    #                            sep='\t', index=False)
    #     with open(args.output_dir + curr_sample_id + '_rna_variants.fasta', 'w') as f:
    #         for curr_element in variant_transcript_sequences:
    #             f.write('>' + curr_element[0] + '\n')
    #             f.write(curr_element[1] + '\n')
    #     logger.info("Finished writing simulated RNA variants files [%i/%i]."
    #                 % (curr_sample_idx + 1, args.num_samples))
    #
    # logger.info("Finished running exacto 'sim-rna-variants' command.")
    #
