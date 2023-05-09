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
and run Exacto 'transcribe' command.
"""


import argparse
from .fasta import *
from .main import *
from .gencode import *
from .default_parameters import *


logger = get_logger(__name__)


def add_cli_simulate_reads_arg_parser(
        sub_parsers
    ) -> argparse._SubParsersAction:
    """
    Adds 'transcribe' parser.

    Parameters
    ----------
    sub_parsers  :   An instance of argparse.ArgumentParser subparsers.

    Returns
    -------
    An instance of argparse.ArgumentParser subparsers.
    """
    parser = sub_parsers.add_parser(
        'transcribe',
        help='Transcribes a genic sequence with DNA variants.'
    )
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        '--reference-genome-fasta-file',
        dest='reference_genome_fasta_file',
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    parser_required.add_argument(
        '--variants-tsv-file',
        dest='variants_tsv_file',
        type=str,
        required=True,
        help="Variants TSV file."
    )
    parser_required.add_argument(
        '--reference-transcripts-gtf-file',
        dest='reference_transcripts_gtf_file',
        type=str,
        required=True,
        help="GENCODE reference transcripts GTF file."
    )
    parser_required.add_argument(
        '--output_transcripts_fasta_file',
        dest='output_transcripts_fasta_file',
        type=str,
        required=True,
        help="Output variants TSV file."
    )


    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser.set_defaults(which='sim-reads')
    return sub_parsers


def run_cli_sim_reads_from_parsed_args(
        args
    ) -> None:
    """
    Run Exacto 'sim-reads' command using parameters from parsed arguments.

    Parameters
    ----------
    args    :   An instance of argparse.ArgumentParser with the following variables:
                fasta_file
                num_gigabases
                output_fastq_file
                read_length_mean
                read_length_stdev
                base_quality_mean
                base_quality_stdev
    """
    logger.info("Started running exacto 'sim-reads' command.")

    # Step 1. Load sequences from the FASTA file
    sequences = read_fasta_file(fasta_file=args.fasta_file)

    # Step 2. Simulate reads
    reads = run_exacto_simulate_reads(
        sequences=sequences,
        output_fastq_gz_file=args.output_fastq_gz_file,
        num_gigabases=args.num_gigabases,
        read_length_mean=args.read_length_mean,
        read_length_stdev=args.read_length_stdev,
        base_quality_mean=args.base_quality_mean,
        base_quality_stdev=args.base_quality_stdev
    )
    #
    # # Step 3. Save to FASTQ file
    # if args.gzip:
    #     logger.info("Started writing simulated reads to FASTQ.GZ file.")
    #     if args.output_fastq_file[-3:] == '.gz':
    #         output_file = args.output_fastq_file
    #     else:
    #         output_file = args.output_fastq_file + '.gz'
    #     with gzip.open(output_file, 'wb') as f:
    #         for read in reads:
    #             f.write(str(read.id + '\n').encode())
    #             f.write(str(read.sequence + '\n').encode())
    #             f.write(str('+' + '\n').encode())
    #             f.write(str(read.base_quality_score_string + '\n').encode())
    #     logger.info("Finished writing simulated reads to FASTQ.GZ file.")
    # else:
    #     logger.info("Started writing simulated reads to FASTQ file.")
    #     with open(args.output_fastq_file, 'w') as f:
    #         for read in reads:
    #             f.write(read.id + '\n')
    #             f.write(read.sequence + '\n')
    #             f.write('+\n')
    #             f.write(read.base_quality_score_string + '\n')
    #     logger.info("Finished writing simulated reads to FASTQ file.")
    #
    logger.info("Finished running exacto 'sim-reads' command.")
