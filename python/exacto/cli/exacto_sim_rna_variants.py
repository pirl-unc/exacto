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


import pysam
from ..constants import *
from ..default_parameters import *
from ..logging import get_logger
from ..main import *
from ..utilities.gencode_utils import *


logger = get_logger(__name__)


def add_exacto_simulate_rna_variants_arg_parser(sub_parsers):
    """
    Adds 'sim-rna-variants' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'sim-rna-variants',
        help='Simulate RNA variants.'
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
        '--gencode_transcripts_gtf_file',
        dest='gencode_transcripts_gtf_file',
        type=str,
        required=True,
        help="GENCODE transcripts GTF file."
    )
    parser_required.add_argument(
        '--transcript_types',
        dest='transcript_types',
        nargs='+',
        required=True,
        help="Transcript types to simulate variants (e.g. 'protein_coding')."
    )
    parser_required.add_argument(
        '--output_variant_transcripts_tsv_file',
        dest='output_variant_transcripts_tsv_file',
        type=str,
        required=True,
        help="Output RNA variants TSV file."
    )
    parser_required.add_argument(
        '--output_variant_transcripts_fasta_file',
        dest='output_variant_transcripts_fasta_file',
        type=str,
        required=True,
        help="Output variant transcripts FASTA file."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--target_regions_tsv_file",
        dest="target_regions_tsv_file",
        type=str,
        required=False,
        help="If this parameter is not specified, all transcripts specified in "
             "--gencode_transcripts_gtf_file will be subject to RNA variant simulation. "
             "If specified, Transcripts within the genomic regions specified in this "
             "file will be simulated for variants. "
             "Expected headers: 'chrom', 'start', 'end'"
    )
    parser_optional.add_argument(
        "--num_snv",
        dest="num_snv",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_SNV,
        required=False,
        help="Number of SNVs to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_SNV
    )
    parser_optional.add_argument(
        "--num_insertion",
        dest="num_insertion",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_INSERTION,
        required=False,
        help="Number of insertions to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_INSERTION
    )
    parser_optional.add_argument(
        "--num_deletion",
        dest="num_deletion",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_DELETION,
        required=False,
        help="Number of deletions to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_DELETION
    )
    parser_optional.add_argument(
        "--num_fusion",
        dest="num_fusion",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_FUSION,
        required=False,
        help="Number of fusions to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_FUSION
    )
    parser_optional.add_argument(
        "--num_inversion",
        dest="num_inversion",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_INVERSION,
        required=False,
        help="Number of inversions to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_INVERSION
    )
    parser_optional.add_argument(
        "--num_intron_retention",
        dest="num_intron_retention",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_INTRON_RETENTION,
        required=False,
        help="Number of intron retentions to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_INTRON_RETENTION
    )
    parser_optional.add_argument(
        "--num_herv",
        dest="num_herv",
        type=int,
        default=SIMULATE_RNA_VARIANTS_NUM_HERV,
        required=False,
        help="Number of human endogenous retroviruses (HERVs) to simulate (default: %i)."
             % SIMULATE_RNA_VARIANTS_NUM_HERV
    )
    parser_optional.add_argument(
        "--herv_regions_tsv_file",
        dest="herv_regions_tsv_file",
        type=str,
        required=False,
        help="This file must be specified if --num_herv is greater than 0. "
             "Expected headers: 'chr', 'start', 'end', 'strand'. "
             "The gEVE Hsap38.txt file can be supplied directly."
    )
    parser_optional.add_argument(
        "--herv_solo_ltr_proportion",
        dest="herv_solo_ltr_proportion",
        type=float,
        default=SIMULATE_RNA_VARIANTS_HERV_PROPORTION_SOLO_LTR,
        required=False,
        help="Proportion of expressed HERVs that only have the solo-LTR sequences expressed (default: %f)."
             % SIMULATE_RNA_VARIANTS_HERV_PROPORTION_SOLO_LTR
    )
    parser_optional.add_argument(
        "--herv_truncated_proportion",
        dest="herv_truncated_proportion",
        type=float,
        default=SIMULATE_RNA_VARIANTS_HERV_PROPORTION_TRUNCATED,
        required=False,
        help="Proportion of expressed HERVs that have truncated sequences (default: %f)."
             % SIMULATE_RNA_VARIANTS_HERV_PROPORTION_TRUNCATED
    )
    parser_optional.add_argument(
        "--herv_chimeric_proportion",
        dest="herv_chimeric_proportion",
        type=float,
        default=SIMULATE_RNA_VARIANTS_HERV_PROPORTION_CHIMERIC,
        required=False,
        help="Proportion of expressed HERVs that are chimeric (neighboring HERVs concatenated; default: %f)."
             % SIMULATE_RNA_VARIANTS_HERV_PROPORTION_CHIMERIC
    )
    parser_optional.add_argument(
        "--herv_chimeric_max_neighboring_distance",
        dest="herv_chimeric_max_neighboring_distance",
        type=int,
        default=SIMULATE_RNA_VARIANTS_HERV_CHIMERIC_MAX_NEIGHBORING_DISTANCE,
        required=False,
        help="Maximum neighboring distance for two HERVs to be simulated as "
             "a chimeric HERV (default: %i)."
             % SIMULATE_RNA_VARIANTS_HERV_CHIMERIC_MAX_NEIGHBORING_DISTANCE
    )
    parser_optional.add_argument(
        "--herv_full_length_proportion",
        dest="herv_full_length_proportion",
        type=float,
        default=SIMULATE_RNA_VARIANTS_HERV_PROPORTION_FULL_LENGTH,
        required=False,
        help="Proportion of expressed HERVs that are full-lengths (default: %f)."
             % SIMULATE_RNA_VARIANTS_HERV_PROPORTION_FULL_LENGTH
    )
    parser_optional.add_argument(
        "--infinite_sites_assumption",
        dest="infinite_sites_assumption",
        type=bool,
        default=SIMULATE_RNA_VARIANTS_INFINITE_SITES_ASSUMPTION,
        required=False,
        help="If true, the simulation enforces infinites sites assumption (default: %r)."
             % SIMULATE_RNA_VARIANTS_INFINITE_SITES_ASSUMPTION
    )
    parser.set_defaults(which='sim-rna-variants')
    return sub_parsers


def run_exacto_simulate_rna_variants_from_parsed_args(args):
    """
    Run Exacto 'sim-rna-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                reference_genome_fasta_file
                gencode_transcripts_gtf_file
                transcript_types
                output_variant_transcripts_tsv_file
                output_variant_transcripts_fasta_file
                target_regions_tsv_file
                num_snv
                num_insertion
                num_deletion
                num_fusion
                num_inversion
                num_intron_retention
                num_herv
                herv_solo_ltr_proportion
                herv_truncated_proportion
                herv_chimeric_proportion
                herv_chimeric_max_neighboring_distance
                herv_full_length_proportion
                infinite_sites_assumption
    """
    # Step 1. Load data
    genome_fasta = pysam.FastaFile(args.reference_genome_fasta_file)
    df_genes, df_transcripts, df_exons = read_gencode_gtf_file(gencode_gtf_file=args.gencode_transcripts_gtf_file)
    if args.target_regions_tsv_file:
        df_target_regions = pd.read_csv(args.target_regions_tsv_file, sep='\t')
    else:
        df_target_regions = None
    if args.num_herv > 0:
        if args.herv_regions_tsv_file:
            df_herv_regions = pd.read_csv(args.herv_regions_tsv_file, sep='\t')
        else:
            logger.error("--herv_regions_tsv_file must be supplied if --num_herv is greater than 0.")
    else:
        df_herv_regions = None

    # Step 2. Simulate RNA variants
    df_transcripts = df_transcripts.loc[df_transcripts['transcript_type'].isin(args.transcript_types),:]
    df_rna_variants, variant_transcript_sequences = run_exacto_simulate_rna_variants(
        genome_fasta=genome_fasta,
        df_genes=df_genes,
        df_transcripts=df_transcripts,
        df_exons=df_exons,
        df_target_regions=df_target_regions,
        df_herv_regions=df_herv_regions,
        num_snv=args.num_snv,
        num_insertion=args.num_insertion,
        num_deletion=args.num_deletion,
        num_fusion=args.num_fusion,
        num_inversion=args.num_inversion,
        num_herv=args.num_herv,
        herv_solo_ltr_proportion=args.herv_solo_ltr_proportion,
        herv_truncated_proportion=args.herv_truncated_proportion,
        herv_chimeric_proportion=args.herv_chimeric_proportion,
        herv_chimeric_max_neighboring_distance=args.herv_chimeric_max_neighboring_distance,
        herv_full_length_proportion=args.herv_full_length_proportion,
        infinite_sites_assumption=args.infinite_sites_assumption
    )


