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
The purpose of this python3 script is to define Exacto default parameters.
"""


# Only keep structural variants with the "precise" tag.
KEEP_ONLY_PRECISE_SV = True

# Only keep structural variants with the following FILTER values.
KEEP_ONLY_FILTER_VALUES = ['PASS']

# Minimum genomic total coverage for a variant position.
MIN_GENOMIC_VARIANT_POSITION_TOTAL_COVERAGE = 7

# Minimum genomic variant reads count.
MIN_GENOMIC_VARIANT_READS_COUNT = 3

# Padding for a genome gapped region.
# The pad is applied to upstream and downstream of a gapped genomic region.
GENOME_GAPPED_REGIONS_PADDING = 100000

# Padding for an excluded structural variant breakpoint.
# The pad is applied to upstream and downstream of the two breakpoints
# of a structural variant to be excluded.
EXCLUDE_SV_PADDING = 20

# Padding for an excluded single-nucleotide variant or
# small insertion / deletion (INDEL).
EXCLUDE_SNV_INDEL_PADDING = 1

# Maximum distance (bases) for merging two structural variants.
MAX_SV_CLUSTER_DISTANCE = 100

# ANNOVAR protocol and corresponding operation
ANNOVAR_PROTOCOL_OPERATION = {
    'refGene': 'g',
    'exac03': 'f',
    '1000g2015aug_eur': 'f',
    '1000g2015aug_eas': 'f',
    '1000g2015aug_sas': 'f',
    'clinvar_20210501': 'f',
    'cosmic96_coding': 'f',
    'avsnp150': 'f',
    'dbnsfp42c': 'f'
}

# Probability of simulating a genic variant
SIMULATE_GENIC_VARIANT_PROBABILITY = 0.2

# Number of single-nucleotide variants to simulate
SIMULATE_NUM_SNV = 100

# Number of small insertions to simulate
SIMULATE_NUM_INSERTION = 100

# Number of small deletions to simulate
SIMULATE_NUM_DELETION = 100

# Enforce infinite sites assumption
SIMULATE_ENFORCE_INFINITE_SITES_ASSUMPTION = True

# Number of cores
NUM_CORES = 4

# Structural variant attributes (union of attributes amongst SV callers)
STRUCTURAL_VARIANT_ATTRIBUTES = {
    'id': 'unknown',                                                           # id
    'variant_calling_method': 'unknown',                                       # variant calling method
    'sequencing_platform': 'unknown',                                          # sequencing platform
    'chr_1': 'unknown',                                                        # chromosome 1
    'pos_1': 'unknown',                                                        # position 1
    'chr_2': 'unknown',                                                        # chromosome 2
    'pos_2': 'unknown',                                                        # position 2
    'ref': 'unknown',                                                          # reference allele
    'alt': 'unknown',                                                          # alternate allele
    'quality_score': 'unknown',                                                # quality score
    'filter': 'unknown',                                                       # filter
    'is_precise': 'unknown',                                                   # is breakpoint precise?
    'sv_type': 'unknown',                                                      # SV type
    'sv_size': 'unknown',                                                      # SV size
    'sv_size_stdev': 'unknown',                                                # SV size standard deviation
    'variant_reads_count': 'unknown',                                          # variant reads count
    'reference_reads_count': 'unknown',                                        # reference reads count
    'total_coverage': 'unknown',                                               # total coverage
    'variant_allele_fraction': 'unknown',                                      # variant allele fraction
    'read_ids': 'unknown',                                                     # read IDs
    'strand': 'unknown',                                                       # strand
    'insertion_sequence': 'unknown',                                           # insertion sequence
    'genotype': 'unknown',                                                     # genotype
    'genotype_quality': 'unknown',                                             # genotype quality
    'sv_pos_stdev': 'unknown',                                                 # SV start position standard deviation
    'coverage': 'unknown',                                                     # coverage (upstream, start, center, end, downstream)
    'query_alignment_length_adjusted_mismatches_mean_count': 'unknown',        # mean number of query alignment length adjusted mismatches of supporting reads
    'support_long': 'unknown',                                                 # number of soft-clipped reads putatively supporting the long insertion SV
    'ci_pos': 'unknown',                                                       # confidence interval around POS for impreicse variants
    'ci_len': 'unknown',                                                       # confidence interval around inserted / deleted material between breakends
    'std_span': 'unknown',                                                     # standard deviation in position of merged SV signatures
    'tandem_duplication_copy_number': 'unknown',                               # copy number of tandem duplication (2 for one additional copy)
    'strand_reads': 'unknown',                                                 # forward and reverse strand reads in each allele
    'repeat_annotation': 'unknown'                                             # repeat annotation
}

# Small variant (SNVs and INDELs) attributes (union of attributes amongst SNV/INDEL callers)
SMALL_VARIANT_ATTRIBUTES = {
    'id': 'unknown',                                                           # id
    'variant_calling_method': 'unknown',                                       # variant calling method
    'sequencing_platform': 'unknown',                                          # sequencing platform
    'chrom': 'unknown',                                                        # chromosome
    'pos': 'unknown',                                                          # position
    'ref': 'unknown',                                                          # reference allele
    'alt': 'unknown',                                                          # alternate allele
    'filter': 'unknown',                                                       # filter
    'quality_score': 'unknown',                                                # quality score
    'variant_type': 'unknown',                                                 # variant type
    'variant_sequence': 'unknown',                                             # variant sequence
    'variant_size': 'unknown',                                                 # variant size
    'tumor_genotype': 'unknown',                                               # tumor genotype
    'tumor_genotype_quality': 'unknown',                                       # tumor genotype quality
    'normal_genotype': 'unknown',                                              # normal genotype
    'normal_genotype_quality': 'unknown',                                      # normal genotype quality
    'tumor_total_coverage': 'unknown',                                         # tumor total coverage
    'tumor_reference_reads_count': 'unknown',                                  # tumor reference reads count
    'tumor_variant_reads_count': 'unknown',                                    # tumor variant reads count
    'normal_total_coverage': 'unknown',                                        # normal total coverage
    'normal_reference_reads_count': 'unknown',                                 # normal reference reads count
    'variant_allele_fraction': 'unknown',                                      # variant allele fraction
    'phred_scale_genotype_likelihoods': 'unknown',                             # phred scale genotype likelihoods
    'allele_specific_strand_bias_table': 'unknown',                            # allele-specific forward/reverse read counts for strand bias tests
    'tumor_strand_bias_fisher_exact_test_component_statistics': 'unknown',     # tumor per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias
    'normal_strand_bias_fisher_exact_test_component_statistics': 'unknown',    # normal per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias
    'haplotype_events': 'unknown',                                             # number of events in this haplotype
    'alt_allele_germline_quality': 'unknown',                                  # phred-scale quality that alt alleles are not germline variants
    'allele_median_base_qualities': 'unknown',                                 # median base quality by allele
    'allele_median_fragment_length': 'unknown',                                # median fragment length by allele
    'allele_median_mapping_quality': 'unknown',                                # median mapping quality by allele
    'median_distance_from_read_end': 'unknown',                                # median distance from end of read
    'negative_log10_odds_artifact': 'unknown',                                 # negative log 10 odds of artifact in normal with same allele fraction as tumor
    'log10_odds_artifact': 'unknown',                                          # normal log 10 odds of artifact in normal with same allele fraction as tumor
    'negative_log_10_population': 'unknown',                                   # negative log 10 population allele frequencies of alt alleles
    'log10_likelihood_ratio_score_variant_exists': 'unknown',                  # log 10 likelihood ratio score of variant existing versus not existing
    'tumor_f1r2_reads_count': 'unknown',                                       # F1R2 pair orientation supporting each allele in tumor
    'tumor_f2r1_reads_count': 'unknown',                                       # F2R1 pair orientation supporting each allele in tumor
    'normal_f1r2_reads_count': 'unknown',                                      # F1R2 pair orientation supporting each allele in normal
    'normal_f2r1_reads_count': 'unknown',                                      # F2R1 pair orientation supporting each allele in normal
    'region_end_position': 'unknown',                                          # End position of the region described int this record
    'non_variant_multisite_block': 'unknown',                                  # on-variant multi-site block. Non-variant blocks are defined independently for each sample. All sites in such a block are constrained to be non-variant, have the same filter value, and have sample values {GQX,DP,DPF} in range [x,y], y <= max(x+3,(x*1.3))
    'snv_contextual_homopolymer_length': 'unknown',                            # SNV contextual homopolymer length
    'cigar': 'unknown',                                                        # CIGAR alignment for each alternate INDEL allele
    'smallest_repeating_sequence_unit': 'unknown',                             # Smallest repeating sequence unit (RU) extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases
    'smallest_repeating_sequence_unit_reference_repeat_count': 'unknown',      # Number of times RU is repeated in reference
    'smallest_repeating_sequence_unit_allele_repeat_count': 'unknown',         # Number of times RU is repeated in indel allele
    'mapping_quality_root_mean_square': 'unknown',                             # Root mean square of mapping quality
    'tumor_genotype_quality_recalibrated': 'unknown',                          # Empirically calibrated genotype quality score for variant sites of tumor, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}
    'normal_genotype_quality_recalibrated': 'unknown',                         # Empirically calibrated genotype quality score for variant sites of normal, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}
    'filtered_basecalls_prior_to_genotyping': 'unknown',                       # Basecalls filtered from input prior to site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block
    'minimum_filtered_basecall_depth': 'unknown',                              # Minimum filtered basecall depth used for site genotyping within a non-variant multi-site block
    'tumor_allelic_depths_forward_strand': 'unknown',                          # Tumor allelic depths on the forward strand
    'tumor_allelic_depths_reverse_strand': 'unknown',                          # Tumor allelic depths on the reverse strand
    'tumor_read_depth_preceding_indel': 'unknown',                             # Read depth associated with INDEL, taken from the site preceding the INDEL in tumor
    'strand_bias_tag': 'unknown'                                               # Strand bias tag
}
