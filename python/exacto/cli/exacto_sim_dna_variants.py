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
and run Exacto 'simulate_somatic_variants' command.
"""


import pysam
from ..constants import *
from ..default_parameters import *
from ..logging import get_logger
from ..main import *


logger = get_logger(__name__)


def add_exacto_simulate_somatic_variants_arg_parser(sub_parsers):
    # """
    # Adds 'simulate_somatic_variants' parser.
    #
    # Parameters
    # ----------
    # sub_parsers  :   An instance of argparse.ArgumentParser subparsers.
    #
    # Returns
    # -------
    # An instance of argparse.ArgumentParser subparsers.
    # """
    # parser = sub_parsers.add_parser(
    #     'sim-somatic-variants',
    #     help='Simulate somatic variants.'
    # )
    # parser._action_groups.pop()
    #
    # # Required arguments
    # parser_required = parser.add_argument_group('required arguments')
    # parser_required.add_argument(
    #     '--gencode_reference_genome_fasta_file',
    #     dest='gencode_reference_genome_fasta_file',
    #     type=str,
    #     required=True,
    #     help="GENCODE reference genome FASTA file."
    # )
    # parser_required.add_argument(
    #     '--gencode_reference_transcripts_gtf_file',
    #     dest='gencode_reference_transcripts_gtf_file',
    #     type=str,
    #     required=True,
    #     help="GENCODE reference transcripts GTF file."
    # )
    # parser_required.add_argument(
    #     '--output_variants_tsv_file',
    #     dest='output_variants_tsv_file',
    #     type=str,
    #     required=True,
    #     help="Output variants TSV file."
    # )
    # parser_required.add_argument(
    #     '--output_augmented_genome_fasta_file',
    #     dest='output_augmented_genome_fasta_file',
    #     type=str,
    #     required=True,
    #     help="Output augmented genome FASTA file."
    # )
    # parser_required.add_argument(
    #     '--output_augmented_transcriptome_fasta_file',
    #     dest='output_augmented_transcriptome_fasta_file',
    #     type=str,
    #     required=True,
    #     help="Output augmented transcriptome FASTA file."
    # )
    #
    # # Optional arguments
    # parser_optional = parser.add_argument_group('optional arguments')
    # parser_optional.add_argument(
    #     "--genic_variant_probability",
    #     dest="genic_variant_probability",
    #     type=float,
    #     default=SIMULATE_GENIC_VARIANT_PROBABILITY,
    #     required=False,
    #     help="Probability of simulating a genic variant (default: %f). "
    #          "The probability of simulating an intergenic variant is "
    #          "therefore 1 - genic_variant_probability (default: %f)."
    #          % (SIMULATE_GENIC_VARIANT_PROBABILITY,
    #             1 - SIMULATE_GENIC_VARIANT_PROBABILITY)
    # )
    # parser_optional.add_argument(
    #     "--num_snv",
    #     dest="num_snv",
    #     type=int,
    #     default=SIMULATE_NUM_SNV,
    #     required=False,
    #     help="Number of SNVs to simulate (default: %i)."
    #          % SIMULATE_NUM_SNV
    # )
    # parser_optional.add_argument(
    #     "--num_insertion",
    #     dest="num_insertion",
    #     type=int,
    #     default=SIMULATE_NUM_INSERTION,
    #     required=False,
    #     help="Number of insertions to simulate (default: %i)."
    #          % SIMULATE_NUM_INSERTION
    # )
    # parser_optional.add_argument(
    #     "--num_deletion",
    #     dest="num_deletion",
    #     type=int,
    #     default=SIMULATE_NUM_DELETION,
    #     required=False,
    #     help="Number of deletions to simulate (default: %i)."
    #          % SIMULATE_NUM_DELETION
    # )
    # parser_optional.add_argument(
    #     "--infinite_sites_assumption",
    #     dest="infinite_sites_assumption",
    #     type=bool,
    #     default=SIMULATE_ENFORCE_INFINITE_SITES_ASSUMPTION,
    #     required=True,
    #     help="Number of deletions to simulate (default: %r)."
    #          % SIMULATE_ENFORCE_INFINITE_SITES_ASSUMPTION
    # )
    #
    # parser.set_defaults(which='simulate_variants')
    # return sub_parsers
    pass


def run_exacto_simulate_variants_from_parsed_args(args):
    """
    Run Exacto 'sim-somatic-variants' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                nucleic_acid_type
                fasta_file
                output_tsv_file
                num_snv
                num_insertion
                num_deletion
    """
    pass
    # fasta = pysam.FastaFile(filename=args.fasta_file)
    # df_variants = run_exacto_simulate_variants(
    #     fasta=fasta
    # )
    # df_variants.to_csv(args.output_tsv_file, sep='\t', index=False)

