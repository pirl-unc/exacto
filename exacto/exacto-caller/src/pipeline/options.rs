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


use exacto_core::prelude::TranslationStrategy;
use exacto_core::prelude::{BaseQuality, MappingQuality, ReadSupport};
use std::collections::HashSet;
use std::string::ToString;


#[derive(Debug, Clone)]
pub struct IdentifyGermlineDNAVariantsOptions {
    pub error_model: SequencingErrorModel,
    pub prior: DNAVariantPriorModel,
    pub calling: DNAVariantCallingOptions,
    pub filtering: DNAVariantFilteringOptions
}

#[derive(Debug, Clone)]
pub struct IdentifySomaticDNAVariantsOptions {
    pub error_model: SequencingErrorModel,
    pub prior: DNAVariantPriorModel,
    pub calling: DNAVariantCallingOptions,
    pub subtraction: DNAControlSubtractionOptions,
    pub filtering: DNAVariantFilteringOptions
}

#[derive(Debug, Clone)]
pub struct IdentifyRNATranscriptVariantsOptions {
    pub calling: RNAVariantCallingOptions,
    pub annotation: RNAAnnotationOptions
}

#[derive(Debug, Clone)]
pub struct SequencingErrorModel {
    /// Probability that a base is a substitution or indel error.
    pub sequencing_error: f64,

    /// Probability of a homopolymer slip per repeat unit.
    pub slippage_prob: f64
}

#[derive(Debug, Clone)]
pub struct DNAVariantPriorModel {
    /// Expected variant allele fraction:
    /// 0.5 for a heterozygous germline site (germline mode).
    /// 0.25 assuming 50% tumor purity (somatic mode).
    pub allele_fraction: f64,

    /// Prior probability that a site is variant, per base.
    pub mutation_rate: f64
}

#[derive(Debug, Clone)]
pub struct DNAVariantCallingOptions {
    /// Region sizes to process in parallel.
    pub chunk_size: u32,

    /// Maximum number of alignments records for a read to be included in variant calling.
    pub max_records: u32,

    /// Minimum mapping quality for a read to be included in variant calling.
    /// If multiple alignments exist for a read name, then the maximum value is compared.
    pub min_mapping_quality: MappingQuality,

    /// Minimum base quality for a base to be included in variant calling.
    pub min_base_quality: BaseQuality,

    /// Maximum insertion normalized edit distance.
    pub max_ins_norm_edit_distance: f64,

    /// Minimum variant size proportion.
    pub min_size_proportion: f64,

    /// Maximum variant records clustering distance.
    pub max_clustering_distance: u32,

    /// Minimum terminal soft-clip insertion length.
    pub min_terminal_soft_clip_ins_len: u32,

    /// Whether to rescue breakpoints from insertions.
    pub bkpt_rescue: bool,

    /// Minimum insertion length to be considered for re-alignment (to rescue breakpoints).
    pub bkpt_rescue_min_ins_len: u32,

    /// Maximum insertion length to be considered for re-alignment (to rescue breakpoints).
    pub bkpt_rescue_max_ins_len: u32,

    /// Search distance for re-alignment (to rescue breakpoints).
    pub bkpt_rescue_search_distance: u32,

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

    /// Minimum span proportion for insertion re-alignment (to rescue breakpoints).
    pub bkpt_rescue_realignment_min_span_proportion: f64,

    /// Partial order alignment match score (to determine consensus insertion sequence).
    pub poa_match_score: i32,

    /// Partial order alignment mismatch score (to determine consensus insertion sequence).
    pub poa_mismatch_score: i32,

    /// Partial order alignment gap open score (to determine consensus insertion sequence).
    pub poa_gap_open_score: i32,

    /// Partial order alignment gap extend score (to determine consensus insertion sequence).
    pub poa_gap_extend_score: i32,

    /// Positions less than this many bases apart have their read depths counted from one query of
    /// the BAM file. Only the number of queries depends on it, not the depths.
    pub read_depth_max_merge_distance: u32
}

#[derive(Debug, Clone)]
pub struct DNAControlSubtractionOptions {
    /// Maximum number of control reads a variant may have and still be called somatic.
    pub max_control_reads: ReadSupport,

    /// Apply infinite sites assumption (ISA).
    pub apply_infinite_sites_assumption: bool
}

#[derive(Debug, Clone)]
pub struct DNAVariantFilteringOptions {
    /// Maximum false positive rate.
    pub max_fpr: f64,

    /// Minimum number of reads for a DNA variant call.
    pub min_reads: usize,

    /// Minimum total depth.
    pub min_total_depth: usize,

    /// Minimum alternate allele fraction.
    pub min_alt_allele_fraction: f64,

    /// Maximum slippage repeat length.
    pub max_slippage_repeat_len: u32,

    /// Minimum homopolymer length.
    pub min_homopolymer_len: u32,

    /// Minimum dinucleotide context length.
    pub min_dinucleotide_context_len: u32
}

#[derive(Debug, Clone)]
pub struct RNAVariantCallingOptions {
    /// Number of reads to process in parallel.
    pub chunk_size: usize,
    
    /// Minimum mapping quality for a read to be included in variant calling.
    /// If multiple alignments exist for a read name, then the maximum value is compared.
    pub min_mapping_quality: MappingQuality,

    /// Minimum terminal soft-clip insertion length.
    pub min_terminal_soft_clip_ins_len: u32,

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

    /// Maximum pieces for insertion re-alignment (to rescue breakpoints).
    /// todo update this
    pub bkpt_rescue_realignment_max_pieces: u32
}

#[derive(Debug, Clone)]
pub struct RNAAnnotationOptions {
    /// Predict nonsense mediated decay.
    pub predict_nonsense_mediated_decay: bool,

    /// Which complete open reading frames to evaluate.
    pub translation_strategy: TranslationStrategy,

    /// Allowed start codons, expressed as uppercase RNA sequences.
    pub start_codons: HashSet<String>,

    /// Nonsense mediated decay distance threshold.
    pub nmd_distance_threshold: u32
}



// =================================================================================================
// Implementations
// =================================================================================================


/// IdentifyGermlineDNAVariantsOptions
impl IdentifyGermlineDNAVariantsOptions {
    pub const DEFAULT: Self = Self {
        error_model: SequencingErrorModel::DEFAULT,
        prior: DNAVariantPriorModel::GERMLINE,
        calling: DNAVariantCallingOptions::DEFAULT,
        filtering: DNAVariantFilteringOptions::GERMLINE
    };

    pub const PACBIO_HIFI: Self = Self {
        error_model: SequencingErrorModel::PACBIO_HIFI,
        filtering: DNAVariantFilteringOptions::GERMLINE_PACBIO_HIFI,
        ..Self::DEFAULT
    };

    pub const ONT: Self = Self {
        error_model: SequencingErrorModel::ONT,
        filtering: DNAVariantFilteringOptions::GERMLINE_ONT,
        ..Self::DEFAULT
    };
}
impl Default for IdentifyGermlineDNAVariantsOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// IdentifySomaticDNAVariantsOptions
impl IdentifySomaticDNAVariantsOptions {
    pub const DEFAULT: Self = Self {
        error_model: SequencingErrorModel::DEFAULT,
        prior: DNAVariantPriorModel::SOMATIC,
        calling: DNAVariantCallingOptions::DEFAULT,
        subtraction: DNAControlSubtractionOptions::DEFAULT,
        filtering: DNAVariantFilteringOptions::SOMATIC
    };

    pub const PACBIO_HIFI: Self = Self {
        error_model: SequencingErrorModel::PACBIO_HIFI,
        filtering: DNAVariantFilteringOptions::SOMATIC_PACBIO_HIFI,
        ..Self::DEFAULT
    };

    pub const ONT: Self = Self {
        error_model: SequencingErrorModel::ONT,
        filtering: DNAVariantFilteringOptions::SOMATIC_ONT,
        ..Self::DEFAULT
    };
}
impl Default for IdentifySomaticDNAVariantsOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// IdentifyRNATranscriptVariantsOptions
impl Default for IdentifyRNATranscriptVariantsOptions {
    fn default() -> Self {
        Self {
            calling: RNAVariantCallingOptions::DEFAULT,
            annotation: RNAAnnotationOptions::default()
        }
    }
}

/// SequencingErrorModel
impl SequencingErrorModel {
    pub const DEFAULT: Self = Self {
        sequencing_error: 0.01,
        slippage_prob: 0.03         // 3x sequencing error
    };

    pub const PACBIO_HIFI: Self = Self {
        sequencing_error: 0.01,
        slippage_prob: 0.03         // 3x sequencing error
    };

    pub const ONT: Self = Self {
        sequencing_error: 0.05,
        slippage_prob: 0.15         // 3x sequencing error
    };
}
impl Default for SequencingErrorModel {
    fn default() -> Self {
        Self::PACBIO_HIFI
    }
}

/// DNAVariantPriorModel
impl DNAVariantPriorModel {
    pub const GERMLINE: Self = Self {
        allele_fraction: 0.5f64,
        mutation_rate: 1e-3f64
    };

    pub const SOMATIC: Self = Self {
        allele_fraction: 0.25f64,   // 50% tumor purity and ploidy of 2
        mutation_rate: 1e-6f64
    };
}

/// DNAVariantCallingOptions
impl DNAVariantCallingOptions {
    pub const DEFAULT: Self = Self {
        chunk_size: 100_000,
        max_records: 7,
        min_mapping_quality: 4,
        min_base_quality: 20,
        max_ins_norm_edit_distance: 0.5f64,
        min_size_proportion: 0.5f64,
        max_clustering_distance: 1_000,
        min_terminal_soft_clip_ins_len: 4,
        bkpt_rescue: true,
        bkpt_rescue_min_ins_len: 100,
        bkpt_rescue_max_ins_len: 500,
        bkpt_rescue_search_distance: 1_000,
        bkpt_rescue_realignment_gap_open_score: -5,
        bkpt_rescue_realignment_gap_extend_score: -1,
        bkpt_rescue_realignment_k: 19,
        bkpt_rescue_realignment_band_width: 100,
        bkpt_rescue_realignment_min_score_fraction: 0.7f64,
        bkpt_rescue_realignment_min_query_coverage: 0.5f64,
        bkpt_rescue_realignment_min_span_proportion: 0.9f64,
        poa_match_score: 0,
        poa_mismatch_score: 4,
        poa_gap_open_score: 6,
        poa_gap_extend_score: 2,
        read_depth_max_merge_distance: 1_000
    };
}
impl Default for DNAVariantCallingOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// DNAControlSubtractionOptions
impl DNAControlSubtractionOptions {
    pub const DEFAULT: Self = Self {
        max_control_reads: 0,
        apply_infinite_sites_assumption: true
    };
}

/// DNAVariantFilteringOptions
impl DNAVariantFilteringOptions {
    pub const GERMLINE: Self = Self {
        max_fpr: 1e-6f64,
        min_reads: 4,
        min_total_depth: 4,
        min_alt_allele_fraction: 0.2f64,
        max_slippage_repeat_len: 30,
        min_homopolymer_len: 4,
        min_dinucleotide_context_len: 4
    };

    pub const SOMATIC: Self = Self {
        max_fpr: 1e-6f64,
        min_reads: 4,
        min_total_depth: 4,
        min_alt_allele_fraction: 0.05f64,
        max_slippage_repeat_len: 30,
        min_homopolymer_len: 4,
        min_dinucleotide_context_len: 4
    };

    pub const GERMLINE_PACBIO_HIFI: Self = Self {
        min_reads: 4,
        min_total_depth: 4,
        ..Self::GERMLINE
    };

    pub const GERMLINE_ONT: Self = Self {
        min_reads: 5,
        min_total_depth: 5,
        ..Self::GERMLINE
    };

    pub const SOMATIC_PACBIO_HIFI: Self = Self {
        min_reads: 4,
        min_total_depth: 4,
        ..Self::SOMATIC
    };

    pub const SOMATIC_ONT: Self = Self {
        min_reads: 5,
        min_total_depth: 5,
        ..Self::SOMATIC
    };
}

/// RNAVariantCallingOptions
impl RNAVariantCallingOptions {
    pub const DEFAULT: Self = Self {
        chunk_size: 10_000,
        min_mapping_quality: 0,
        min_terminal_soft_clip_ins_len: 4,
        bkpt_rescue: true,
        bkpt_rescue_min_ins_len: 100,
        bkpt_rescue_realignment_gap_open_score: -5,
        bkpt_rescue_realignment_gap_extend_score: -1,
        bkpt_rescue_realignment_k: 11,
        bkpt_rescue_realignment_band_width: 50,
        bkpt_rescue_realignment_min_score_fraction: 0.8f64,
        bkpt_rescue_realignment_min_query_coverage: 0.1f64,
        bkpt_rescue_realignment_min_placed_fraction: 0.9f64,
        bkpt_rescue_realignment_max_pieces: 3
    };
}
impl Default for RNAVariantCallingOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// RNAAnnotationOptions
impl Default for RNAAnnotationOptions {
    fn default() -> Self {
        Self {
            predict_nonsense_mediated_decay: true,
            translation_strategy: TranslationStrategy::LongestORF,
            start_codons: HashSet::from(["AUG".to_string()]),
            nmd_distance_threshold: 50
        }
    }
}