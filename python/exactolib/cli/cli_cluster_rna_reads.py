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
and run Exacto 'cluster-rna-reads' command.
"""


import argparse
import os

from ..main import *
from ..utilities import *


logger = get_logger(__name__)


def add_cli_cluster_rna_reads_arg_parser(sub_parsers) -> argparse._SubParsersAction:
    """
    Add 'cluster-rna-reads' parser.

    Parameters:
        sub_parsers     :  argparse.ArgumentParser subparsers.

    Returns:
        sub_parsers     :   argparse.ArgumentParser subparsers
    """
    parser = sub_parsers.add_parser('cluster-rna-reads', help='Cluster RNA reads (full-length long RNA reads).')
    parser._action_groups.pop()

    # Required arguments
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument(
        "--bam-file",
        dest="bam_file",
        type=str,
        required=True,
        help="Input BAM file."
    )
    parser_required.add_argument(
        "--bai-file",
        dest="bai_file",
        type=str,
        required=True,
        help="Input BAM.BAI file."
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
    parser_optional.add_argument(
        "--preset",
        dest="preset",
        type=str,
        choices=list(CLUSTER_RNA_READS_PRESETS.keys()),
        default=None,
        required=False,
        help="Error-profile preset that fills typical defaults "
             "for --expected-sequencing-error, --expected-slippage-probability, "
             "--min-reads-per-cluster, --min-reads, and --min-total-depth. "
             "Choices: 'pb' (PacBio HiFi), "
             "'ont' (Oxford Nanopore), or 'corrected' (exacto correct-rna-reads output, "
             "for pass-2 re-clustering). Any explicit parameter wins over the preset."
    )
    parser_optional.add_argument(
        "--analyte-type",
        dest="analyte_type",
        type=str,
        choices=['cdna', 'rna'],
        default=CLUSTER_RNA_READS_ANALYTE_TYPE,
        required=False,
        help="Analyte sequenced: 'cdna' (made by reverse transcription) or 'rna' (direct RNA "
             "sequencing, which turns reverse transcriptase template-switch filtering off) "
             "(default: %s)." % CLUSTER_RNA_READS_ANALYTE_TYPE
    )
    parser_optional.add_argument(
        "--remove-unspliced-rnas",
        dest="remove_unspliced_rnas",
        type=str2bool,
        default=CLUSTER_RNA_READS_REMOVE_UNSPLICED_RNAS,
        required=False,
        help="Remove unspliced RNAs, as remove-unspliced-rnas does, before clustering; pass 'no' "
             "when the input BAM is already its output (default: %s)."
             % CLUSTER_RNA_READS_REMOVE_UNSPLICED_RNAS
    )
    parser_optional.add_argument(
        "--max-batch-reads",
        dest="max_batch_reads",
        type=int,
        default=CLUSTER_RNA_READS_MAX_BATCH_READS,
        required=False,
        help="Maximum number of reads to process in a batch (default: %s)."
             % CLUSTER_RNA_READS_MAX_BATCH_READS
    )
    parser_optional.add_argument(
        "--max-locus-gap",
        dest="max_locus_gap",
        type=int,
        default=CLUSTER_RNA_READS_MAX_LOCUS_GAP,
        required=False,
        help="Maximum locus gap for initial read clustering (default: %s)."
             % CLUSTER_RNA_READS_MAX_LOCUS_GAP
    )
    parser_optional.add_argument(
        "--unspliced-bin-size",
        dest="unspliced_bin_size",
        type=int,
        default=CLUSTER_RNA_READS_UNSPLICED_BIN_SIZE,
        required=False,
        help="Unspliced read bin size for initial read clustering (default: %s)."
             % CLUSTER_RNA_READS_UNSPLICED_BIN_SIZE
    )
    parser_optional.add_argument(
        "--min-mapping-quality",
        dest="min_mapping_quality",
        type=int,
        default=CLUSTER_RNA_READS_MIN_MAPPING_QUALITY,
        required=False,
        help="Minimum mapping quality (default: %i)."
             % CLUSTER_RNA_READS_MIN_MAPPING_QUALITY
    )
    parser_optional.add_argument(
        "--min-terminal-softclip-length",
        dest="min_terminal_softclip_length",
        type=int,
        default=CLUSTER_RNA_READS_MIN_TERMINAL_SOFTCLIP_LENGTH,
        required=False,
        help="Terminal soft-clip runs shorter than this many bases emit no variant record; "
             "1-2 bp terminal clips are aligner anchor jitter, not junction or insertion "
             "evidence. Set 0 to disable (default: %i)."
             % CLUSTER_RNA_READS_MIN_TERMINAL_SOFTCLIP_LENGTH
    )
    parser_optional.add_argument(
        "--min-reads-per-cluster",
        dest="min_reads_per_cluster",
        type=int,
        default=None,
        required=False,
        help="Minimum number of reads per splice-junction cluster. This is the junction-cluster "
             "read floor only; variant calls are gated separately by --min-reads and "
             "--min-total-depth (default: preset-dependent, pb 3 / ont 4 / "
             "corrected 3; %i with no preset). A flat floor, and the only read-support test a "
             "splice-junction cluster has to pass."
             % CLUSTER_RNA_READS_MIN_READS_PER_CLUSTER
    )
    parser_optional.add_argument(
        "--max-records",
        dest="max_records",
        type=int,
        default=CLUSTER_RNA_READS_MAX_RECORDS,
        required=False,
        help="Maximum number of records. Read names having more than this value will be excluded (default: %i)."
             % CLUSTER_RNA_READS_MAX_RECORDS
    )
    parser_optional.add_argument(
        "--expected-sequencing-error",
        dest="expected_sequencing_error",
        type=float,
        default=None,
        required=False,
        help="Expected sequencing error rate (default: %f with no preset, %f with --preset pb, %f with --preset ont)."
             % (CLUSTER_RNA_READS_EXPECTED_SEQUENCING_ERROR,
                CLUSTER_RNA_READS_PRESETS["pb"]["expected_sequencing_error"],
                CLUSTER_RNA_READS_PRESETS["ont"]["expected_sequencing_error"])
    )
    parser_optional.add_argument(
        "--expected-slippage-probability",
        dest="expected_slippage_probability",
        type=float,
        default=None,
        required=False,
        help="Expected slippage probability (default: %f with no preset, %f with --preset pb, %f with --preset ont). "
             "Suggested: 2x the expected sequencing error rate."
             % (CLUSTER_RNA_READS_EXPECTED_SLIPPAGE_PROBABILITY,
                CLUSTER_RNA_READS_PRESETS["pb"]["expected_slippage_probability"],
                CLUSTER_RNA_READS_PRESETS["ont"]["expected_slippage_probability"])
    )
    parser_optional.add_argument(
        "--max-fpr",
        dest="max_fpr",
        type=float,
        default=CLUSTER_RNA_READS_MAX_FPR,
        required=False,
        help="Maximum false positive rate of a variant call within a cluster. A call needs enough "
             "supporting reads that P(X >= reads | site depth, --expected-sequencing-error) falls "
             "below this; inside a repeat the error rate is the slippage probability "
             "(--expected-slippage-probability, --max-slippage-repeat-length). Retained introns "
             "and novel splice junctions, including unannotated consecutive-junction pairs, "
             "use the same beta-binomial error test at local read "
             "depth, with mean and dispersion both set to --expected-sequencing-error. "
             "The separate --min-reads-per-cluster floor also applies (default: %s)."
             % CLUSTER_RNA_READS_MAX_FPR
    )
    parser_optional.add_argument(
        "--max-slippage-repeat-length",
        dest="max_slippage_repeat_length",
        type=int,
        default=CLUSTER_RNA_READS_MAX_SLIPPAGE_REPEAT_LENGTH,
        required=False,
        help="Maximum slippage repeat length (default: %i)."
             % CLUSTER_RNA_READS_MAX_SLIPPAGE_REPEAT_LENGTH
    )
    parser_optional.add_argument(
        "--min-size-proportion",
        dest="min_size_proportion",
        type=float,
        default=CLUSTER_RNA_READS_MIN_SIZE_PROPORTION,
        required=False,
        help="Minimum size proportion between two insertions for them to cluster (default: %f). "
             "Size proportion = smaller variant size / longer variant size."
             % CLUSTER_RNA_READS_MIN_SIZE_PROPORTION
    )
    parser_optional.add_argument(
        "--max-ins-norm-edit-distance",
        dest="max_ins_norm_edit_distance",
        type=float,
        default=CLUSTER_RNA_READS_MAX_INS_NORM_EDIT_DISTANCE,
        required=False,
        help="Maximum normalized edit distance between two insertion sequences for them to "
             "cluster (default: %f)."
             % CLUSTER_RNA_READS_MAX_INS_NORM_EDIT_DISTANCE
    )
    parser_optional.add_argument(
        "--max-intrachromosomal-distance",
        dest="max_intrachromosomal_distance",
        type=int,
        default=CLUSTER_RNA_READS_MAX_INTRACHROMOSOMAL_DISTANCE,
        required=False,
        help="Maximum distance between two insertions for them to cluster (default: %i)."
             % CLUSTER_RNA_READS_MAX_INTRACHROMOSOMAL_DISTANCE
    )
    parser_optional.add_argument(
        "--dna-variants-tsv-file",
        dest="dna_variants_tsv_file",
        type=str,
        required=False,
        help="DNA variants TSV file (an output of exacto call-germline-dna-vars or "
             "call-somatic-dna-vars). When supplied, an RNA variant call that matches a listed DNA "
             "variant exactly (chromosomes, positions, operations and sequence) passes without the "
             "depth, read-support and template-switch tests. Cluster membership and "
             "--min-reads-per-cluster are unaffected."
    )
    parser_optional.add_argument(
        "--allowed-variants-tsv-file",
        dest="allowed_variants_tsv_file",
        type=str,
        required=False,
        help="Allowed variants TSV file: the _exacto_rna_clusters_variants_passed.tsv of a "
             "pass-1 exacto cluster-rna-reads run, for pass 2 of the two-pass RNA pipeline. When "
             "supplied, a variant call that matches no listed variant (same type and chromosomes, "
             "both positions within --allowed-variant-max-distance, an insertion's sequence within "
             "--max-ins-norm-edit-distance) fails as 'not_in_allowed_list' and is not phased on. "
             "A listed call still faces every other test."
    )
    parser_optional.add_argument(
        "--allowed-variant-max-distance",
        dest="allowed_variant_max_distance",
        type=int,
        default=CLUSTER_RNA_READS_ALLOWED_VARIANT_MAX_DISTANCE,
        required=False,
        help="Largest distance in bases between a call's positions and those of the allowed "
             "variant it matches; correction and realignment can move a variant between passes "
             "(default: %i)."
             % CLUSTER_RNA_READS_ALLOWED_VARIANT_MAX_DISTANCE
    )
    parser_optional.add_argument(
        "--mec-max-k",
        dest="mec_max_k",
        type=int,
        default=CLUSTER_RNA_READS_MEC_MAX_K,
        required=False,
        help="Maximum number of clusters to identify at the single isoform level for minimum error correction (default: %i)."
             % CLUSTER_RNA_READS_MEC_MAX_K
    )
    parser_optional.add_argument(
        "--mec-num-restarts",
        dest="mec_num_restarts",
        type=int,
        default=CLUSTER_RNA_READS_MEC_NUM_RESTARTS,
        required=False,
        help="Number of restarts for the minimum error correction algorithm (default: %i)."
             % CLUSTER_RNA_READS_MEC_NUM_RESTARTS
    )
    parser_optional.add_argument(
        "--mec-max-iter",
        dest="mec_max_iter",
        type=int,
        default=CLUSTER_RNA_READS_MEC_MAX_ITER,
        required=False,
        help="Maximum number of iterations the minimum error correction algorithm (default: %i)."
             % CLUSTER_RNA_READS_MEC_MAX_ITER
    )
    parser_optional.add_argument(
        "--mec-seed",
        dest="mec_seed",
        type=int,
        default=CLUSTER_RNA_READS_MEC_SEED,
        required=False,
        help="Seed for the minimum error correction algorithm (default: %i)."
             % CLUSTER_RNA_READS_MEC_SEED
    )
    parser_optional.add_argument(
        "--ts-min-homology",
        dest="ts_min_homology",
        type=int,
        default=CLUSTER_RNA_READS_TS_MIN_HOMOLOGY,
        required=False,
        help="Junction-spanning flank homology (bases) that flags a junction on its own, "
             "unless it carries the GT-AG spliceosomal signature (default: %i)."
             % CLUSTER_RNA_READS_TS_MIN_HOMOLOGY
    )
    parser_optional.add_argument(
        "--ts-soft-min-homology",
        dest="ts_soft_min_homology",
        type=int,
        default=CLUSTER_RNA_READS_TS_SOFT_MIN_HOMOLOGY,
        required=False,
        help="Homology at or above this flags only with corroboration: breakpoint "
             "dispersion, fold-back geometry, or hub promiscuity (default: %i)."
             % CLUSTER_RNA_READS_TS_SOFT_MIN_HOMOLOGY
    )
    parser_optional.add_argument(
        "--ts-max-breakpoint-dispersion",
        dest="ts_max_breakpoint_dispersion",
        type=int,
        default=CLUSTER_RNA_READS_TS_MAX_BREAKPOINT_DISPERSION,
        required=False,
        help="Pooled member breakpoints spreading more than this many bases beyond the "
             "homology interval count as dispersion corroboration (default: %i)."
             % CLUSTER_RNA_READS_TS_MAX_BREAKPOINT_DISPERSION
    )
    parser_optional.add_argument(
        "--num-threads",
        dest="num_threads",
        type=int,
        default=CLUSTER_RNA_READS_NUM_THREADS,
        required=False,
        help="Number of threads (default: %i)."
             % CLUSTER_RNA_READS_NUM_THREADS
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
        "--max-unspliced-locus-gap",
        dest="max_unspliced_locus_gap",
        type=int,
        default=CLUSTER_RNA_READS_MAX_UNSPLICED_LOCUS_GAP,
        required=False,
        help="Maximum allowed distance between two unspliced RNA reads for them to share a locus "
             "(default: %i)."
             % CLUSTER_RNA_READS_MAX_UNSPLICED_LOCUS_GAP
    )
    parser_optional.add_argument(
        "--soft-clip-removal",
        dest="soft_clip_removal",
        type=str2bool,
        default=CLUSTER_RNA_READS_SOFT_CLIP_REMOVAL,
        required=False,
        help="Whether to remove the insertions made of a terminal soft clip that spells the "
             "reference across an intron. Such a clip is the end of a degraded read reaching into "
             "the neighbouring exon, too short for the aligner to splice. Also resolve a clip "
             "that spells a fusion, breakpoint or translocation junction other reads align, as "
             "those reads spell it past the junction (see --soft-clip-min-partner-bases), as that junction, "
             "and move its read to the one cluster carrying it (default: %s)."
             % CLUSTER_RNA_READS_SOFT_CLIP_REMOVAL
    )
    parser_optional.add_argument(
        "--soft-clip-max-boundary-distance",
        dest="soft_clip_max_boundary_distance",
        type=int,
        default=CLUSTER_RNA_READS_SOFT_CLIP_MAX_BOUNDARY_DISTANCE,
        required=False,
        help="Maximum distance in bases between the end of the alignment and the boundary of the "
             "intron or junction a terminal soft clip is compared across (for --soft-clip-removal) "
             "(default: %i)."
             % CLUSTER_RNA_READS_SOFT_CLIP_MAX_BOUNDARY_DISTANCE
    )
    parser_optional.add_argument(
        "--soft-clip-bases-per-edit",
        dest="soft_clip_bases_per_edit",
        type=int,
        default=CLUSTER_RNA_READS_SOFT_CLIP_BASES_PER_EDIT,
        required=False,
        help="Number of bases of a terminal soft clip for every edit allowed between the clip and "
             "the reference across the intron, or the supporting reads across the junction (for --soft-clip-removal). Set 0 to allow no edit "
             "(default: %i)."
             % CLUSTER_RNA_READS_SOFT_CLIP_BASES_PER_EDIT
    )
    parser_optional.add_argument(
        "--soft-clip-min-partner-bases",
        dest="soft_clip_min_partner_bases",
        type=int,
        default=CLUSTER_RNA_READS_SOFT_CLIP_MIN_PARTNER_BASES,
        required=False,
        help="Minimum number of bases of a terminal soft clip, past the untemplated bases of a "
             "fusion, breakpoint or translocation junction, that must spell the junction's partner "
             "arm for the clip to be resolved as that junction (for --soft-clip-removal). A shorter "
             "stretch can match an unrelated arm by chance (default: %i)."
             % CLUSTER_RNA_READS_SOFT_CLIP_MIN_PARTNER_BASES
    )
    parser_optional.add_argument(
        "--unplaced-tail-removal",
        dest="unplaced_tail_removal",
        type=str2bool,
        default=CLUSTER_RNA_READS_UNPLACED_TAIL_REMOVAL,
        required=False,
        help="Whether to remove the stretch of a read the aligner left unplaced around a short "
             "aligned block: a large insertion, a splice junction, the block and a terminal soft "
             "clip, where the block is shorter than the insertion and the clip. The junction into "
             "the block and the variant records of the stretch are removed (default: %s)."
             % CLUSTER_RNA_READS_UNPLACED_TAIL_REMOVAL
    )
    parser_optional.add_argument(
        "--unplaced-tail-min-ins-len",
        dest="unplaced_tail_min_ins_len",
        type=int,
        default=CLUSTER_RNA_READS_UNPLACED_TAIL_MIN_INS_LEN,
        required=False,
        help="Minimum insertion length for a stretch to be removed as an unplaced tail (for "
             "--unplaced-tail-removal) (default: %i)."
             % CLUSTER_RNA_READS_UNPLACED_TAIL_MIN_INS_LEN
    )
    parser_optional.add_argument(
        "--chunk-size",
        dest="chunk_size",
        type=int,
        default=CLUSTER_RNA_READS_CHUNK_SIZE,
        required=False,
        help="Number of reads to process in parallel per chunk (default: %i)."
             % CLUSTER_RNA_READS_CHUNK_SIZE
    )
    parser_optional.add_argument(
        "--bkpt-rescue",
        dest="bkpt_rescue",
        type=str2bool,
        default=CLUSTER_RNA_READS_BKPT_RESCUE,
        required=False,
        help="Whether to rescue breakpoints from insertions by re-aligning long insertion sequences "
             "(default: %s)."
             % CLUSTER_RNA_READS_BKPT_RESCUE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-min-ins-len",
        dest="bkpt_rescue_min_ins_len",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_MIN_INS_LEN,
        required=False,
        help="Minimum insertion length to be considered for re-alignment (to rescue breakpoints) "
             "(default: %i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_MIN_INS_LEN
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-open-score",
        dest="bkpt_rescue_realignment_gap_open_score",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE,
        required=False,
        help="Gap open score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_OPEN_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-gap-extend-score",
        dest="bkpt_rescue_realignment_gap_extend_score",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE,
        required=False,
        help="Gap extend score for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_GAP_EXTEND_SCORE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-k",
        dest="bkpt_rescue_realignment_k",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_K,
        required=False,
        help="K-mer size for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_K
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-band-width",
        dest="bkpt_rescue_realignment_band_width",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH,
        required=False,
        help="Band width for insertion re-alignment (to rescue breakpoints) (default: %i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_BAND_WIDTH
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-score-fraction",
        dest="bkpt_rescue_realignment_min_score_fraction",
        type=float,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION,
        required=False,
        help="Minimum score fraction for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_SCORE_FRACTION
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-query-coverage",
        dest="bkpt_rescue_realignment_min_query_coverage",
        type=float,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE,
        required=False,
        help="Minimum query coverage for insertion re-alignment (to rescue breakpoints) (default: %f)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_QUERY_COVERAGE
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-min-placed-fraction",
        dest="bkpt_rescue_realignment_min_placed_fraction",
        type=float,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION,
        required=False,
        help="Minimum placed fraction for insertion re-alignment (to rescue breakpoints) (default: "
             "%f)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MIN_PLACED_FRACTION
    )
    parser_optional.add_argument(
        "--bkpt-rescue-realignment-max-pieces",
        dest="bkpt_rescue_realignment_max_pieces",
        type=int,
        default=CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES,
        required=False,
        help="Maximum number of pieces for insertion re-alignment (to rescue breakpoints) (default: "
             "%i)."
             % CLUSTER_RNA_READS_BKPT_RESCUE_REALIGNMENT_MAX_PIECES
    )
    parser_optional.add_argument(
        "--poa-match-score",
        dest="poa_match_score",
        type=int,
        default=CLUSTER_RNA_READS_POA_MATCH_SCORE,
        required=False,
        help="Partial order alignment match score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CLUSTER_RNA_READS_POA_MATCH_SCORE
    )
    parser_optional.add_argument(
        "--poa-mismatch-score",
        dest="poa_mismatch_score",
        type=int,
        default=CLUSTER_RNA_READS_POA_MISMATCH_SCORE,
        required=False,
        help="Partial order alignment mismatch score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CLUSTER_RNA_READS_POA_MISMATCH_SCORE
    )
    parser_optional.add_argument(
        "--poa-gap-open-score",
        dest="poa_gap_open_score",
        type=int,
        default=CLUSTER_RNA_READS_POA_GAP_OPEN_SCORE,
        required=False,
        help="Partial order alignment gap open score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CLUSTER_RNA_READS_POA_GAP_OPEN_SCORE
    )
    parser_optional.add_argument(
        "--poa-gap-extend-score",
        dest="poa_gap_extend_score",
        type=int,
        default=CLUSTER_RNA_READS_POA_GAP_EXTEND_SCORE,
        required=False,
        help="Partial order alignment gap extend score (to determine consensus insertion sequence) "
             "(default: %i)."
             % CLUSTER_RNA_READS_POA_GAP_EXTEND_SCORE
    )
    parser_optional.add_argument(
        "--min-reads",
        dest="min_reads",
        type=int,
        default=None,
        required=False,
        help="Minimum number of reads supporting a variant call (default: preset-dependent, "
             "pb 3 / ont 4 / corrected 3; %i with no preset)."
             % CLUSTER_RNA_READS_MIN_READS
    )
    parser_optional.add_argument(
        "--min-total-depth",
        dest="min_total_depth",
        type=int,
        default=None,
        required=False,
        help="Minimum total depth at a variant site (default: preset-dependent, "
             "pb 3 / ont 4 / corrected 3; %i with no preset)."
             % CLUSTER_RNA_READS_MIN_TOTAL_DEPTH
    )
    parser_optional.add_argument(
        "--min-homopolymer-len",
        dest="min_homopolymer_len",
        type=int,
        default=CLUSTER_RNA_READS_MIN_HOMOPOLYMER_LEN,
        required=False,
        help="Minimum homopolymer length for a variant to be treated as homopolymer-context (default: "
             "%i)."
             % CLUSTER_RNA_READS_MIN_HOMOPOLYMER_LEN
    )
    parser_optional.add_argument(
        "--min-dinucleotide-context-len",
        dest="min_dinucleotide_context_len",
        type=int,
        default=CLUSTER_RNA_READS_MIN_DINUCLEOTIDE_CONTEXT_LEN,
        required=False,
        help="Minimum dinucleotide repeat context length for a variant to be treated as "
             "dinucleotide-context (default: %i)."
             % CLUSTER_RNA_READS_MIN_DINUCLEOTIDE_CONTEXT_LEN
    )
    parser_optional.add_argument(
        "--template-switch-flank",
        dest="template_switch_flank",
        type=int,
        default=CLUSTER_RNA_READS_TEMPLATE_SWITCH_FLANK,
        required=False,
        help="Flank length in bases examined on each side of a junction for reverse transcriptase "
             "template-switch homology (default: %i)."
             % CLUSTER_RNA_READS_TEMPLATE_SWITCH_FLANK
    )
    parser_optional.add_argument(
        "--template-switch-foldback-min-stem",
        dest="template_switch_foldback_min_stem",
        type=int,
        default=CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MIN_STEM,
        required=False,
        help="Minimum stem length for the reverse transcriptase template-switch fold-back signature "
             "(default: %i)."
             % CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MIN_STEM
    )
    parser_optional.add_argument(
        "--template-switch-foldback-max-loop-len",
        dest="template_switch_foldback_max_loop_len",
        type=int,
        default=CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_LOOP_LEN,
        required=False,
        help="Maximum loop length for the reverse transcriptase template-switch fold-back signature "
             "(default: %i)."
             % CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_LOOP_LEN
    )
    parser_optional.add_argument(
        "--template-switch-foldback-max-distance",
        dest="template_switch_foldback_max_distance",
        type=int,
        default=CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_DISTANCE,
        required=False,
        help="Maximum distance (stem + loop + stem) for the reverse transcriptase template-switch "
             "fold-back signature (default: %i)."
             % CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_MAX_DISTANCE
    )
    parser_optional.add_argument(
        "--template-switch-foldback-slack",
        dest="template_switch_foldback_slack",
        type=int,
        default=CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_SLACK,
        required=False,
        help="Slack in bases for the reverse transcriptase template-switch fold-back signature "
             "(default: %i)."
             % CLUSTER_RNA_READS_TEMPLATE_SWITCH_FOLDBACK_SLACK
    )

    parser.set_defaults(which='cluster-rna-reads')
    return sub_parsers


def run_cli_cluster_rna_reads_vars_from_parsed_args(args) -> None:
    """
    Run Exacto 'cluster-rna-reads' command using parameters from parsed arguments.
    """
    preset_values: dict[str, int | float] = CLUSTER_RNA_READS_PRESETS.get(args.preset, {})
    expected_sequencing_error = (
        args.expected_sequencing_error
        if args.expected_sequencing_error is not None
        else preset_values.get("expected_sequencing_error", CLUSTER_RNA_READS_EXPECTED_SEQUENCING_ERROR)
    )
    expected_slippage_probability = (
        args.expected_slippage_probability
        if args.expected_slippage_probability is not None
        else preset_values.get("expected_slippage_probability", CLUSTER_RNA_READS_EXPECTED_SLIPPAGE_PROBABILITY)
    )
    min_reads_per_cluster = (
        args.min_reads_per_cluster
        if args.min_reads_per_cluster is not None
        else preset_values.get("min_reads_per_cluster", CLUSTER_RNA_READS_MIN_READS_PER_CLUSTER)
    )
    min_reads = (
        args.min_reads
        if args.min_reads is not None
        else preset_values.get("min_reads", CLUSTER_RNA_READS_MIN_READS)
    )
    min_total_depth = (
        args.min_total_depth
        if args.min_total_depth is not None
        else preset_values.get("min_total_depth", CLUSTER_RNA_READS_MIN_TOTAL_DEPTH)
    )

    if args.dna_variants_tsv_file is None:
        dna_variants_tsv_file = ""
    else:
        dna_variants_tsv_file = args.dna_variants_tsv_file

    os.makedirs(args.output_dir, exist_ok=True)

    cluster_rna_reads(
        bam_file=args.bam_file,
        bai_file=args.bai_file,
        reference_genome_fasta_file=args.reference_genome_fasta_file,
        reference_gene_annotation_file=args.reference_gene_annotation_file,
        reference_gene_annotation_assembly=args.reference_gene_annotation_assembly,
        reference_gene_annotation_source=args.reference_gene_annotation_source,
        reference_gene_annotation_version=args.reference_gene_annotation_version,
        output_dir=args.output_dir,
        output_prefix=args.output_prefix,
        analyte_type=args.analyte_type,
        remove_unspliced_rnas=args.remove_unspliced_rnas,
        max_batch_reads=args.max_batch_reads,
        max_locus_gap=args.max_locus_gap,
        unspliced_bin_size=args.unspliced_bin_size,
        min_mapping_quality=args.min_mapping_quality,
        min_terminal_softclip_length=args.min_terminal_softclip_length,
        min_reads_per_cluster=min_reads_per_cluster,
        max_records=args.max_records,
        expected_sequencing_error=expected_sequencing_error,
        expected_slippage_probability=expected_slippage_probability,
        max_fpr=args.max_fpr,
        max_slippage_repeat_length=args.max_slippage_repeat_length,
        min_size_proportion=args.min_size_proportion,
        max_ins_norm_edit_distance=args.max_ins_norm_edit_distance,
        max_intrachromosomal_distance=args.max_intrachromosomal_distance,
        dna_variants_tsv_file=dna_variants_tsv_file,
        allowed_variants_tsv_file="" if args.allowed_variants_tsv_file is None else args.allowed_variants_tsv_file,
        allowed_variant_max_distance=args.allowed_variant_max_distance,
        mec_max_k=args.mec_max_k,
        mec_num_restarts=args.mec_num_restarts,
        mec_max_iter=args.mec_max_iter,
        mec_seed=args.mec_seed,
        ts_min_homology=args.ts_min_homology,
        ts_soft_min_homology=args.ts_soft_min_homology,
        ts_max_breakpoint_dispersion=args.ts_max_breakpoint_dispersion,
        max_unspliced_locus_gap=args.max_unspliced_locus_gap,
        chunk_size=args.chunk_size,
        soft_clip_removal=args.soft_clip_removal,
        soft_clip_max_boundary_distance=args.soft_clip_max_boundary_distance,
        soft_clip_bases_per_edit=args.soft_clip_bases_per_edit,
        soft_clip_min_partner_bases=args.soft_clip_min_partner_bases,
        unplaced_tail_removal=args.unplaced_tail_removal,
        unplaced_tail_min_ins_len=args.unplaced_tail_min_ins_len,
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
        poa_match_score=args.poa_match_score,
        poa_mismatch_score=args.poa_mismatch_score,
        poa_gap_open_score=args.poa_gap_open_score,
        poa_gap_extend_score=args.poa_gap_extend_score,
        min_reads=min_reads,
        min_total_depth=min_total_depth,
        min_homopolymer_len=args.min_homopolymer_len,
        min_dinucleotide_context_len=args.min_dinucleotide_context_len,
        template_switch_flank=args.template_switch_flank,
        template_switch_foldback_min_stem=args.template_switch_foldback_min_stem,
        template_switch_foldback_max_loop_len=args.template_switch_foldback_max_loop_len,
        template_switch_foldback_max_distance=args.template_switch_foldback_max_distance,
        template_switch_foldback_slack=args.template_switch_foldback_slack,
        num_threads=args.num_threads,
        temp_dir=args.temp_dir
    )
