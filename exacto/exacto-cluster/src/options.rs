// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use exacto_core::prelude::AnalyteType;


#[derive(Debug, Clone)]
pub struct ClusterRNAReadsOptions {
    /// If AnalyteType is AnalyteType::RNA, then reverse transcriptase template-switch
    /// artifact filtering is turned OFF (default: AnalyteType::CDNA).
    pub analyte_type: AnalyteType,
    pub seed: usize,
    pub remove_unspliced_rnas: bool,
    pub max_batch_reads: usize,

    /// Two junctions of one read are one locus when they sit on one chromosome closer than this.
    /// A read whose junctions fall in several loci is a chimera or a fusion; its loci are never
    /// unioned with one another, so a handful of such reads cannot merge the expressed genome.
    pub max_locus_gap: u32,

    /// A junction-less read overlapping no annotated gene is keyed by the bin holding its start.
    /// Bins are never linked to one another: linking through straddling reads chains every bin of
    /// a chromosome together at library depth.
    pub unspliced_bin_size: u32,

    pub error_model: exacto_caller::prelude::SequencingErrorModel,
    pub junction: JunctionClusteringOptions,
    pub calling: RNAVariantCallingOptions,
    pub filtering: RNAVariantFilteringOptions,
    pub phasing: PhasingOptions
}

#[derive(Debug, Clone)]
pub struct JunctionClusteringOptions {
    /// Minimum number of reads for a splice junction cluster.
    pub min_reads: usize,

    /// Maximum allowed distance between two unspliced RNA reads.
    pub max_unspliced_locus_gap: u32
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RNAVariantCallingOptions {
    /// Number of reads to process in parallel.
    pub chunk_size: usize,

    /// Maximum number of alignments records for a read to be included in variant calling.
    pub max_records: usize,

    /// Minimum mapping quality for a read to be included in variant calling.
    /// todo update this: If multiple alignments exist for a read name, then...?
    pub min_mapping_quality: u16,

    /// Maximum variant records clustering distance.
    pub max_clustering_distance: u32,

    /// Maximum insertion normalized edit distance.
    pub max_ins_norm_edit_distance: f64,

    /// Minimum variant size proportion.
    pub min_size_proportion: f64,

    /// Minimum terminal soft-clipped insertion length.
    pub min_terminal_soft_clip_ins_len: u32,

    /// Whether to remove terminal clips spelling the reference across an intron (the end of a
    /// degraded read reaching into the neighbouring exon), and to resolve those spelling a fusion,
    /// breakpoint or translocation junction other reads align as that junction.
    pub soft_clip_removal: bool,

    /// Maximum distance between the end of the alignment and the intron or junction boundary a
    /// terminal soft clip is compared across.
    pub soft_clip_max_boundary_distance: u32,

    /// Number of bases of a terminal soft clip for every edit allowed between the clip and the
    /// reference across the intron, or the bases a supporting read spells across the junction.
    /// 0 allows no edit.
    pub soft_clip_bases_per_edit: u32,

    /// Fewest bases of a terminal soft clip, past the untemplated bases of a fusion, breakpoint
    /// or translocation junction, that must spell the junction's partner arm for the clip to be
    /// resolved as that junction. A shorter stretch can match an unrelated arm by chance.
    pub soft_clip_min_partner_bases: u32,

    /// Whether to remove the stretch of a read the aligner left unplaced around a short aligned
    /// block: a large insertion, a splice junction, the block and a terminal soft clip.
    pub unplaced_tail_removal: bool,

    /// Minimum insertion length for a stretch to be removed as an unplaced tail.
    pub unplaced_tail_min_ins_len: u32,

    /// Whether to rescue breakpoints from insertions.
    pub bkpt_rescue: bool,

    /// Minimum insertion length to be considered for re-alignment (to rescue breakpoints).
    pub bkpt_rescue_min_ins_len: u32,

    /// Gap open score for insertion re-alignment (to rescue breakpoints).
    pub bkpt_rescue_realignment_gap_open_score: i32,

    /// Gap extend score for insertion re-alignment (to rescue breakpoints).
    pub bkpt_rescue_realignment_gap_extend_score: i32,

    /// K for insertion re-alignment (to rescue breakpoints).
    pub bkpt_rescue_realignment_k: u32,

    /// Band width for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_band_width: u32,

    /// Minimum score fraction for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_min_score_fraction: f64,

    /// Minimum score fraction for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_min_query_coverage: f64,

    /// Minimum placed fraction for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_min_placed_fraction: f64,

    /// Maximum number of pieces for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_max_pieces: u32,

    /// Partial order alignment match score (to determine consensus insertion sequence).
    pub poa_match_score: i32,

    /// Partial order alignment mismatch score (to determine consensus insertion sequence).
    pub poa_mismatch_score: i32,

    /// Partial order alignment gap open score (to determine consensus insertion sequence).
    pub poa_gap_open_score: i32,

    /// Partial order alignment gap extend score (to determine consensus insertion sequence).
    pub poa_gap_extend_score: i32
}

#[derive(Debug, Clone)]
pub struct RNAVariantFilteringOptions {
    /// Maximum false positive rate.
    pub max_fpr: f64,

    /// Minimum number of reads for an RNA variant call.
    pub min_reads: usize,

    /// Minimum total depth (within a splice junction cluster).
    pub min_total_depth: usize,

    /// Maximum slippage repeat length.
    pub max_slippage_repeat_len: u32,

    /// Minimum homopolymer length.
    pub min_homopolymer_len: u32,

    /// Minimum dinucleotide context length.
    pub min_dinucleotide_context_len: u32,

    /// Reverse transcriptase template switch flank.
    /// todo: update this doc.
    pub template_switch_flank: u32,

    /// Minimum homology for reverse transcriptase template switch junction.
    /// todo: update this doc.
    pub template_switch_junction_min_homology: u32,

    /// Minimum soft homology for reverse transcriptase template switch junction.
    /// todo: update this doc.
    pub template_switch_junction_soft_min_homology: u32,

    /// Maximum breakpoint dispersion for reverse transcriptase template switch junction.
    /// todo: update this doc.
    pub template_switch_junction_max_breakpoint_dispersion: u32,

    /// Minimum stem for reverse transcriptase template switch foldback signature.
    /// todo: update this doc.
    pub template_switch_foldback_min_stem: u32,

    /// Maximum loop length for reverse transcriptase template switch foldback signature.
    /// todo: update this doc.
    pub template_switch_foldback_max_loop_len: u32,

    /// Maximum distance (stem + loop + stem) for reverse transcriptase template switch foldback signature.
    /// todo: update this doc.
    pub template_switch_foldback_max_distance: u32,

    /// Slack for reverse transcriptase template switch foldback signature.
    /// todo: update this doc.
    pub template_switch_foldback_slack: u32,

    /// Largest distance, in bases, between a call's positions and those of the allowed variant
    /// it matches, when the run is restricted to an allow list. A variant's positions can move
    /// between two passes: correction and realignment re-spell it.
    pub allowed_variant_max_distance: u32
}

#[derive(Debug, Clone)]
pub struct PhasingOptions {
    /// Minimum error correction maximum k.
    pub mec_max_k: usize,

    /// Minimum error correction number of restarts.
    pub mec_num_restart: usize,

    /// Minimum error correction number of maximum iterations.
    pub mec_max_iter: usize
}


/// ClusterRNAReadsOptions
impl ClusterRNAReadsOptions {
    pub const DEFAULT: Self = Self {
        seed: 42,
        analyte_type: AnalyteType::CDNA,
        remove_unspliced_rnas: true,
        max_batch_reads: 100_000,
        max_locus_gap: 1_000_000,
        unspliced_bin_size: 100_000,
        error_model: exacto_caller::prelude::SequencingErrorModel::DEFAULT,
        junction: JunctionClusteringOptions::DEFAULT,
        calling: RNAVariantCallingOptions::DEFAULT,
        filtering: RNAVariantFilteringOptions::DEFAULT,
        phasing: PhasingOptions::DEFAULT
    };

    pub const PACBIO_HIFI: Self = Self {
        error_model: exacto_caller::prelude::SequencingErrorModel::PACBIO_HIFI,
        junction: JunctionClusteringOptions::PACBIO_HIFI,
        filtering: RNAVariantFilteringOptions::PACBIO_HIFI,
        ..Self::DEFAULT
    };

    pub const ONT: Self = Self {
        error_model: exacto_caller::prelude::SequencingErrorModel::ONT,
        junction: JunctionClusteringOptions::ONT,
        filtering: RNAVariantFilteringOptions::ONT,
        ..Self::DEFAULT
    };
}
impl Default for ClusterRNAReadsOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// JunctionClusteringOptions
impl JunctionClusteringOptions {
    pub const DEFAULT: Self = Self {
        min_reads: 3,
        max_unspliced_locus_gap: 0
    };

    pub const PACBIO_HIFI: Self = Self {
        min_reads: 3,
        ..Self::DEFAULT
    };

    pub const ONT: Self = Self {
        min_reads: 4,
        ..Self::DEFAULT
    };
}
impl Default for JunctionClusteringOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// VariantCallingOptions
impl RNAVariantCallingOptions {
    pub const DEFAULT: Self = Self {
        chunk_size: 10_000,
        max_records: 7,
        min_mapping_quality: 0,
        max_clustering_distance: 1_000,
        max_ins_norm_edit_distance: 0.5f64,
        min_size_proportion: 0.5f64,
        min_terminal_soft_clip_ins_len: 4,
        soft_clip_removal: true,
        soft_clip_max_boundary_distance: 3,
        soft_clip_bases_per_edit: 8,
        soft_clip_min_partner_bases: 20,
        unplaced_tail_removal: true,
        unplaced_tail_min_ins_len: 100,
        bkpt_rescue: true,
        bkpt_rescue_min_ins_len: 100,
        bkpt_rescue_realignment_gap_open_score: -5,
        bkpt_rescue_realignment_gap_extend_score: -1,
        bkpt_rescue_realignment_k: 11,
        bkpt_rescue_realignment_band_width: 50,
        bkpt_rescue_realignment_min_score_fraction: 0.8f64,
        bkpt_rescue_realignment_min_query_coverage: 0.1f64,
        bkpt_rescue_realignment_min_placed_fraction: 0.9f64,
        bkpt_rescue_realignment_max_pieces: 3,
        poa_match_score: 0,
        poa_mismatch_score: 4,
        poa_gap_open_score: 6,
        poa_gap_extend_score: 2
    };
}
impl Default for RNAVariantCallingOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// VariantFilteringOptions
impl RNAVariantFilteringOptions {
    pub const DEFAULT: Self = Self {
        max_fpr: 1e-6f64,
        min_reads: 3,
        min_total_depth: 3,
        max_slippage_repeat_len: 30,
        min_homopolymer_len: 4,
        min_dinucleotide_context_len: 4,
        template_switch_flank: 20,
        template_switch_junction_min_homology: 8,
        template_switch_junction_soft_min_homology: 5,
        template_switch_junction_max_breakpoint_dispersion: 2,
        template_switch_foldback_min_stem: 5,
        template_switch_foldback_max_loop_len: 20,
        template_switch_foldback_max_distance: 1_000,
        template_switch_foldback_slack: 20,
        allowed_variant_max_distance: 10
    };

    pub const PACBIO_HIFI: Self = Self {
        min_reads: 3,
        min_total_depth: 3,
        ..Self::DEFAULT
    };

    pub const ONT: Self = Self {
        min_reads: 4,
        min_total_depth: 4,
        ..Self::DEFAULT
    };
}
impl Default for RNAVariantFilteringOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// PhasingOptions
impl PhasingOptions {
    pub const DEFAULT: Self = Self {
        mec_max_k: 7,
        mec_num_restart: 1_000,
        mec_max_iter: 10_000
    };
}
impl Default for PhasingOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}
