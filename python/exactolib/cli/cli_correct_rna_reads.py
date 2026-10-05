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
and run Exacto 'correct-rna-reads' command.
"""


import argparse
from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_correct_rna_reads_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'correct-rna-reads' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser(
        'correct-rna-reads',
        help="Error-correct RNA reads against the cluster each was assigned to.")
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file. Must be the same alignment exacto cluster-rna-reads was run on. Reads are recovered "
             "from their primary records in their original orientation, with their base quality scores."
    )
    parser_required.add_argument(
        "--clusters-tsv-file",
        dest="clusters_tsv_file",
        type=str,
        required=True,
        help="Clusters TSV file (an output of exacto cluster-rna-reads). "
             "Expected columns: 'cluster_id', 'read_name'."
    )
    parser_required.add_argument(
        "--cluster-reference-transcripts-tsv-file",
        dest="cluster_reference_transcripts_tsv_file",
        type=str,
        required=True,
        help="Cluster reference transcripts TSV file (an output of exacto cluster-rna-reads)."
    )
    parser_required.add_argument(
        "--cluster-splice-junctions-tsv-file",
        dest="cluster_splice_junctions_tsv_file",
        type=str,
        required=True,
        help="Cluster splice junctions TSV file (an output of exacto cluster-rna-reads)."
    )
    parser_required.add_argument(
        "--cluster-variants-tsv-file",
        dest="cluster_variants_tsv_file",
        type=str,
        required=True,
        help="Cluster variants TSV file (an output of exacto cluster-rna-reads). Each read is rewritten "
             "to the reference plus its cluster's called alleles: a base, insertion or deletion that "
             "disagrees with them is corrected, and at a called variant the read takes the called "
             "allele, so this file is what separates real variation from sequencing error."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-file",
        dest="reference_gene_annotation_file",
        type=str,
        default="",
        required=True,
        help="Reference gene annotation file."
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
        required=True,
        help="Output prefix."
    )

    # Optional arguments
    parser_optional = parser.add_argument_group('optional arguments')
    parser_optional.add_argument(
        "--max-correctable-event-len",
        dest="max_correctable_event_len",
        type=int,
        default=CORRECT_RNA_READS_MAX_CORRECTABLE_EVENT_LEN,
        required=False,
        help="Maximum length of an uncalled insertion or deletion that correction will rewrite; longer "
             "ones are left as sequenced. Called alleles are written at any length, and terminal soft "
             "clips are kept unless their cluster called them as an insertion (default: %i)."
             % CORRECT_RNA_READS_MAX_CORRECTABLE_EVENT_LEN
    )
    parser_optional.add_argument(
        "--corrected-base-quality",
        dest="corrected_base_quality",
        type=int,
        default=CORRECT_RNA_READS_CORRECTED_BASE_QUALITY,
        required=False,
        help="Base quality assigned to a corrected base (default: %i)."
             % CORRECT_RNA_READS_CORRECTED_BASE_QUALITY
    )
    parser_optional.add_argument(
        "--trim-transcript-ends",
        dest="trim_transcript_ends",
        action="store_true",
        default=CORRECT_RNA_READS_TRIM_TRANSCRIPT_ENDS,
        required=False,
        help="Cut a read end that runs past an annotated exon into the intron back to the exon, when "
             "that exon is not the first or last exon of the reference transcripts that carry the "
             "read's terminal splice junction. A called variant in the removed span keeps the end. "
             "Such an end is also what cluster-rna-reads counts as a retained intron, so trimming "
             "removes retained introns from the corrected reads. Requires "
             "--reference-gene-annotation-file."
    )
    parser_optional.add_argument(
        "--reference-gene-annotation-source",
        dest="reference_gene_annotation_source",
        type=str,
        default=str(GeneAnnotationSource.GENCODE),
        required=False,
        help="Reference gene annotation source (choices: %s)." %
             ','.join([str(GeneAnnotationSource.GENCODE)])
    )
    parser_optional.add_argument(
        "--reference-gene-annotation-assembly",
        dest="reference_gene_annotation_assembly",
        type=str,
        default="",
        required=False,
        help="Reference gene annotation assembly (e.g. 'hg38')."
    )
    parser_optional.add_argument(
        "--reference-gene-annotation-version",
        dest="reference_gene_annotation_version",
        type=str,
        default="",
        required=False,
        help="Reference gene annotation version (e.g. 'v41')."
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CORRECT_RNA_READS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CORRECT_RNA_READS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--chunk-size",
        dest="chunk_size",
        type=int,
        default=CORRECT_RNA_READS_CHUNK_SIZE,
        required=False,
        help="Number of reads corrected per chunk; a cluster is never split across chunks "
             "(default: %i)." % CORRECT_RNA_READS_CHUNK_SIZE
    )

    parser.set_defaults(which='correct-rna-reads')
    return sub_parsers


def run_cli_correct_rna_reads_from_parsed_args(args) -> None:
    """
    Run Exacto 'correct-rna-reads' command using parameters from parsed arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    correct_rna_reads(
        bam_file=args.bam_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_source=GeneAnnotationSource(args.reference_gene_annotation_source),
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        clusters_tsv_file=args.clusters_tsv_file,
        cluster_reference_transcripts_tsv_file=args.cluster_reference_transcripts_tsv_file,
        cluster_splice_junctions_tsv_file=args.cluster_splice_junctions_tsv_file,
        cluster_variants_tsv_file=args.cluster_variants_tsv_file,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        max_correctable_event_len=args.max_correctable_event_len,
        corrected_base_quality=args.corrected_base_quality,
        trim_transcript_ends=args.trim_transcript_ends,
        num_threads=args.num_threads,
        chunk_size=args.chunk_size
    )
