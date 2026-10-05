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
and run Exacto 'call-rna-transcript-vars' command.
"""


import argparse

from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_call_rna_transcript_vars_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'call-rna-transcript-vars' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('call-rna-transcript-vars', help='Call RNA variants in a long-read RNA-seq assembled transcriptome BAM file.')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file of assembled transcripts."
    )
    parser_required.add_argument(
        "--reference-genome-fasta-file",
        dest="reference_genome_fasta_file",
        type=str,
        required=True,
        help="Reference genome FASTA file."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-file",
        dest="reference_gene_annotation_file",
        type=str,
        required=True,
        help="Reference gene annotation file."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-source",
        dest="reference_gene_annotation_source",
        type=str,
        required=True,
        help="Reference gene annotation source (choices: %s)." %
             ','.join([str(GeneAnnotationSource.GENCODE)])
    )
    parser_required.add_argument(
        "--reference-gene-annotation-assembly",
        dest="reference_gene_annotation_assembly",
        type=str,
        required=True,
        help="Reference gene annotation assembly (e.g. 'hg38')."
    )
    parser_required.add_argument(
        "--reference-gene-annotation-version",
        dest="reference_gene_annotation_version",
        type=str,
        required=True,
        help="Reference gene annotation version (e.g. 'v41')."
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
    # parser_optional.add_argument(
    #     "--cluster-variants-tsv-file",
    #     dest="cluster_variants_tsv_file",
    #     type=str,
    #     default="",
    #     required=False,
    #     help="Passed cluster variants TSV file (*_exacto_rna_clusters_variants_passed.tsv). Expected columns: 'cluster_id', 'chromosome_1', 'position_1', 'strand_1', 'operation_type_1', 'chromosome_2', 'position_2', 'strand_2', 'operation_type_2', 'sequence', 'variant_type'. "
    #          "This file is one of the outputs from exacto cluster-rna-reads. "
    #          "When supplied, each transcript keeps only the variants its own cluster called, so every read name in --bam-file must be a cluster ID "
    #          "(default: keep every variant)."
    # )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CALL_RNA_VARS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CALL_RNA_VARS_NUM_THREADS
    )
    parser_optional.add_argument(
        "--chunk-size",
        dest="chunk_size",
        type=int,
        default=CALL_RNA_VARS_CHUNK_SIZE,
        required=False,
        help="Number of reads processed per chunk (default: %i)." % CALL_RNA_VARS_CHUNK_SIZE
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=CALL_RNA_VARS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % CALL_RNA_VARS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-terminal-softclip-length",
        dest="min_terminal_softclip_length",
        type=int,
        default=CALL_RNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        required=False,
        help="Terminal soft-clip runs shorter than this many bases emit no variant record; "
             "1-2 bp terminal clips are aligner anchor jitter, not junction or insertion "
             "evidence. Set 0 to disable (default: %i)."
             % CALL_RNA_VARS_MIN_TERMINAL_SOFTCLIP_LENGTH
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
        "--bkpt-rescue",
        dest="bkpt_rescue",
        type=str2bool,
        default=CALL_RNA_VARS_BKPT_RESCUE,
        required=False,
        help="Whether to rescue breakpoints from insertions by re-aligning long insertion sequences "
             "(default: %s)."
             % CALL_RNA_VARS_BKPT_RESCUE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-min-ins-len",
        dest="bkpt_rescue_min_ins_len",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_MIN_INS_LEN,
        required=False,
        help="Minimum insertion length to be considered for re-alignment (to rescue breakpoints) "
             "(default: %i)."
             % CALL_RNA_VARS_BKPT_RESCUE_MIN_INS_LEN
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-open-score",
        dest="bkpt_rescue_realignment_gap_open_score",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        required=False,
        help="Gap open score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-extend-score",
        dest="bkpt_rescue_realignment_gap_extend_score",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        required=False,
        help="Gap extend score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-k",
        dest="bkpt_rescue_realignment_k",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_K,
        required=False,
        help="K-mer size for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_K
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-band-width",
        dest="bkpt_rescue_realignment_band_width",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        required=False,
        help="Band width for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-score-fraction",
        dest="bkpt_rescue_realignment_min_score_fraction",
        type=float,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        required=False,
        help="Minimum score fraction for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-query-coverage",
        dest="bkpt_rescue_realignment_min_query_coverage",
        type=float,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        required=False,
        help="Minimum query coverage for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-placed-fraction",
        dest="bkpt_rescue_realignment_min_placed_fraction",
        type=float,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION,
        required=False,
        help="Minimum placed fraction for insertion re-alignment (to rescue breakpoints) (default: "
             "%f)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-max-pieces",
        dest="bkpt_rescue_realignment_max_pieces",
        type=int,
        default=CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES,
        required=False,
        help="Maximum number of pieces for insertion re-alignment (to rescue breakpoints) (default: "
             "%i)."
             % CALL_RNA_VARS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES
    )
    parser_optional.add_argument(
        "--predict-nonsense-mediated-decay",
        dest="predict_nonsense_mediated_decay",
        type=str2bool,
        default=CALL_RNA_VARS_PREDICT_NONSENSE_MEDIATED_DECAY,
        required=False,
        help="Whether to predict nonsense-mediated decay for the open reading frames of each transcript "
             "(default: %s)."
             % CALL_RNA_VARS_PREDICT_NONSENSE_MEDIATED_DECAY
    )
    parser_optional.add_argument(
        "--translation-strategy",
        dest="translation_strategy",
        type=str,
        choices=['longest_orf', 'all_orfs'],
        default=CALL_RNA_VARS_TRANSLATION_STRATEGY,
        required=False,
        help="Translation strategy: which complete open reading frames of each transcript to evaluate "
             "when predicting nonsense-mediated decay (default: %s)."
             % CALL_RNA_VARS_TRANSLATION_STRATEGY
    )
    parser_optional.add_argument(
        "--start-codons",
        dest="start_codons",
        type=str,
        nargs='+',
        default=CALL_RNA_VARS_START_CODONS,
        required=False,
        help="One or more start codons. Pass multiple as: --start-codons AUG GUG CUG (default: %s)."
             % ' '.join(CALL_RNA_VARS_START_CODONS)
    )
    parser_optional.add_argument(
        "--nmd-distance-threshold",
        dest="nmd_distance_threshold",
        type=int,
        default=CALL_RNA_VARS_NMD_DISTANCE_THRESHOLD,
        required=False,
        help="Distance in nucleotides upstream of the last exon-exon junction beyond which a stop "
             "codon marks the transcript as NMD-sensitive (default: %i)."
             % CALL_RNA_VARS_NMD_DISTANCE_THRESHOLD
    )
    parser_optional.add_argument(
        "--dna-variants-tsv-files",
        dest="dna_variants_tsv_files",
        type=str,
        nargs='+',
        default=[],
        required=False,
        help="One or more DNA variant TSV files written by 'call-germline-dna-vars' or "
             "'call-somatic-dna-vars' (column 'origin' required). An RNA variant with the same "
             "chromosomes, operations, variant type and sequence as a DNA variant, each position "
             "within --dna-variant-match-buffer bases (strand not compared), carries the DNA "
             "variant's origin in the 'origin' column of the RNA variants TSV file (default: none)."
    )
    parser_optional.add_argument(
        "--dna-variant-match-buffer",
        dest="dna_variant_match_buffer",
        type=int,
        default=CALL_RNA_VARS_DNA_VARIANT_MATCH_BUFFER,
        required=False,
        help="Bases each position of a matching DNA variant may differ by (default: %i)."
             % CALL_RNA_VARS_DNA_VARIANT_MATCH_BUFFER
    )
    parser.set_defaults(which='call-rna-transcript-vars')
    return sub_parsers


def run_cli_call_rna_transcript_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'call-rna-transcript-vars' command using parameters from parsed arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    identify_rna_transcript_variants(
        bam_file=args.bam_file,
        reference_genome_fasta_file=args.reference_genome_fasta_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_source=GeneAnnotationSource(args.reference_gene_annotation_source),
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        # cluster_variants_tsv_file=args.cluster_variants_tsv_file,
        min_mapping_quality=args.min_mapping_quality,
        min_terminal_softclip_length=args.min_terminal_softclip_length,
        bkpt_rescue=args.bkpt_rescue,
        bkpt_rescue_min_ins_len=args.bkpt_rescue_min_ins_len,
        bkpt_rescue_realignment_gap_open_score=args.bkpt_rescue_realignment_gap_open_score,
        bkpt_rescue_realignment_gap_extend_score=args.bkpt_rescue_realignment_gap_extend_score,
        bkpt_rescue_realignment_k=args.bkpt_rescue_realignment_k,
        bkpt_rescue_realignment_band_width=args.bkpt_rescue_realignment_band_width,
        bkpt_rescue_realignment_min_score_fraction=args.bkpt_rescue_realignment_min_score_fraction,
        bkpt_rescue_realignment_min_query_coverage=args.bkpt_rescue_realignment_min_query_coverage,
        bkpt_rescue_realignment_min_placed_fraction=args.bkpt_rescue_realignment_min_placed_fraction,
        bkpt_rescue_realignment_max_pieces=args.bkpt_rescue_realignment_max_pieces,
        predict_nonsense_mediated_decay=args.predict_nonsense_mediated_decay,
        translation_strategy=args.translation_strategy,
        start_codons=args.start_codons,
        nmd_distance_threshold=args.nmd_distance_threshold,
        dna_variants_tsv_files=args.dna_variants_tsv_files,
        dna_variant_match_buffer=args.dna_variant_match_buffer,
        num_threads=args.num_threads,
        chunk_size=args.chunk_size,
        temp_dir=args.temp_dir
    )
