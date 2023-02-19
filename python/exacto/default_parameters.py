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


"""convert"""
# Structural variant attributes (union of attributes amongst SV callers)
STRUCTURAL_VARIANT_ATTRIBUTES = {
    'sample_id': 'unknown',                                                    # sample ID
    'variant_id': 'unknown',                                                   # variant ID
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
    'total_depth': 'unknown',                                                  # total depth
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
    'repeat_annotation': 'unknown'                                             # repeat annotations
}

# Small variant (SNVs and INDELs) attributes (union of attributes amongst SNV/INDEL callers)
SMALL_VARIANT_ATTRIBUTES = {
    'sample_id': '',                                                           # sample ID
    'tumor_sample_id': '',                                                     # tumor sample ID
    'normal_sample_id': '',                                                    # normal sample ID
    'variant_id': 'unknown',                                                   # variant ID
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
    'genotype': 'unknown',                                                     # genotype
    'genotype_quality': 'unknown',                                             # genotype quality
    'genotype_quality_recalibrated': 'unknown',                                # empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}
    'filtered_basecalls_prior_to_genotyping': 'unknown',                       # basecalls filtered from input prior to site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block
    'normal_genotype': 'unknown',                                              # normal genotype
    'normal_genotype_quality': 'unknown',                                      # normal genotype quality
    'normal_genotype_quality_recalibrated': 'unknown',                         # empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}
    'normal_filtered_basecalls_prior_to_genotyping': 'unknown',                # basecalls filtered from input prior to site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block
    'total_depth': 'unknown',                                                  # total depth
    'reference_reads_count': 'unknown',                                        # reference reads count
    'variant_reads_count': 'unknown',                                          # variant reads count
    'variant_allele_fraction': 'unknown',                                      # variant allele fraction
    'normal_total_depth': 'unknown',                                           # normal total depth
    'normal_reference_reads_count': 'unknown',                                 # normal reference reads count
    'phred_scale_genotype_likelihoods': 'unknown',                             # phred-scale genotype likelihoods
    'normal_phred_scale_genotype_likelihoods': 'unknown',                      # normal phred-scale genotype likelihoods
    'allele_specific_strand_bias_table': 'unknown',                            # allele-specific forward/reverse read counts for strand bias tests
    'strand_bias_fisher_exact_test_component_statistics': 'unknown',           # per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias
    'normal_strand_bias_fisher_exact_test_component_statistics': 'unknown',    # normal per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias
    'haplotype_events': 'unknown',                                             # number of events in this haplotype
    'alt_allele_germline_quality': 'unknown',                                  # phred-scale quality that alt alleles are not germline variants
    'allele_median_base_qualities': 'unknown',                                 # median base quality by allele
    'allele_median_fragment_length': 'unknown',                                # median fragment length by allele
    'allele_median_mapping_quality': 'unknown',                                # median mapping quality by allele
    'median_distance_from_read_end': 'unknown',                                # median distance from end of read
    'negative_log10_odds_artifact': 'unknown',                                 # negative log 10 odds of artifact in normal with same allele fraction
    'log10_odds_artifact': 'unknown',                                          # normal log 10 odds of artifact in normal with same allele fraction
    'negative_log_10_population': 'unknown',                                   # negative log 10 population allele frequencies of alt alleles
    'log10_likelihood_ratio_score_variant_exists': 'unknown',                  # log10 likelihood ratio score of variant existing versus not existing
    'f1r2_reads_count': 'unknown',                                             # F1R2 pair orientation supporting each allele
    'f2r1_reads_count': 'unknown',                                             # F2R1 pair orientation supporting each allele
    'normal_f1r2_reads_count': 'unknown',                                      # F1R2 pair orientation supporting each allele in normal
    'normal_f2r1_reads_count': 'unknown',                                      # F2R1 pair orientation supporting each allele in normal
    'region_end_position': 'unknown',                                          # end position of the region described int this record
    'snv_contextual_homopolymer_length': 'unknown',                            # SNV contextual homopolymer length
    'cigar': 'unknown',                                                        # CIGAR alignment for each alternate INDEL allele
    'smallest_repeating_sequence_unit': 'unknown',                             # Smallest repeating sequence unit (RU) extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases
    'smallest_repeating_sequence_unit_reference_repeat_count': 'unknown',      # number of times RU is repeated in reference
    'smallest_repeating_sequence_unit_allele_repeat_count': 'unknown',         # number of times RU is repeated in indel allele
    'mapping_quality_root_mean_square': 'unknown',                             # root-mean-square of mapping quality
    'normal_genotype_quality_recalibrated': 'unknown',                         # empirically calibrated genotype quality score for variant sites of normal, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}
    'minimum_filtered_basecall_depth': 'unknown',                              # minimum filtered basecall depth used for site genotyping within a non-variant multi-site block
    'allelic_depths_forward_strand': 'unknown',                                # allelic depths on the forward strand
    'allelic_depths_reverse_strand': 'unknown',                                # allelic depths on the reverse strand
    'normal_allelic_depths_forward_strand': 'unknown',                         # normal allelic depths on the forward strand
    'normal_allelic_depths_reverse_strand': 'unknown',                         # normal allelic depths on the reverse strand
    'read_depth_preceding_indel': 'unknown',                                   # read depth associated with INDEL, taken from the site preceding the INDEL
    'strand_bias': 'unknown',                                                  # strand bias
    'normal_strand_bias': 'unknown'                                            # normal strand bias
}


"""refine"""
# Only keep structural variants with the "precise" tag.
KEEP_ONLY_PRECISE_SV = True

# Only keep structural variants with the following FILTER values.
KEEP_ONLY_FILTER_VALUES = ['PASS']

# Minimum genomic total depth for a variant position.
MIN_GENOMIC_VARIANT_POSITION_TOTAL_DEPTH = 7

# Minimum genomic variant reads count.
MIN_GENOMIC_VARIANT_READS_COUNT = 3

# Padding for a genome gapped region.
# The pad is applied to upstream and downstream of a gapped genomic region.
GENOME_GAPPED_REGIONS_PADDING = 100000

# Padding for an excluded structural variant breakpoint.
# The pad is applied to upstream and downstream of the two breakpoints
# of a structural variant to be excluded.
EXCLUDE_SV_PADDING = 20

# Enforce variant type checking.
ENFORCE_VARIANT_TYPE_CHECK = True

# Number of processes
NUM_PROCESSES_REFINE = 4


"""merge"""
# Maximum distance (bases) for merging two structural variants.
MAX_SV_CLUSTER_DISTANCE = 10

# Maximum distance (bases) for merging two small variants (i.e. insertions and deletions).
MAX_SMALL_VARIANT_CLUSTER_DISTANCE = 1


"""annotate"""
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


"""sim-dna-variants"""
# Probability of simulating a genic variant
SIMULATE_GENIC_VARIANT_PROBABILITY = 0.2


"""sim-rna-variants"""
# Number of single-nucleotide variants to simulate
SIMULATE_RNA_VARIANTS_NUM_SNV = 300

# Number of small insertions to simulate
SIMULATE_RNA_VARIANTS_NUM_INSERTION = 50

# Insertion size mean
SIMULATE_RNA_VARIANTS_INSERTION_SIZE_MEAN = 100

# Insertion size standard deviation
SIMULATE_RNA_VARIANTS_INSERTION_SIZE_STDEV = 50

# Number of small deletions to simulate
SIMULATE_RNA_VARIANTS_NUM_DELETION = 50

# Deletion size mean
SIMULATE_RNA_VARIANTS_DELETION_MEAN = 100

# Deletion size standard deviation
SIMULATE_RNA_VARIANTS_DELETION_STDEV = 50

# Number of fusion genes to simulate
SIMULATE_RNA_VARIANTS_NUM_FUSION = 5

# Number of inversions to simulate
SIMULATE_RNA_VARIANTS_NUM_INVERSION = 5

# Number of intron retentions to simulate
SIMULATE_RNA_VARIANTS_NUM_INTRON_RETENTION = 5

# Number of HERVs to simulate
SIMULATE_RNA_VARIANTS_NUM_HERV = 20

# Proportion of solo-LTR HERV
SIMULATE_RNA_VARIANTS_HERV_PROPORTION_SOLO_LTR = 0.586 # She et al., Genome Biology 2022

# Proportion of truncated HERV
SIMULATE_RNA_VARIANTS_HERV_PROPORTION_TRUNCATED = 0.23 # She et al., Genome Biology 2022

# Proportion of chimeric HERV
SIMULATE_RNA_VARIANTS_HERV_PROPORTION_CHIMERIC = 0.132 # She et al., Genome Biology 2022

# Chimeric HERV maximum neighboring distance
SIMULATE_RNA_VARIANTS_HERV_CHIMERIC_MAX_NEIGHBORING_DISTANCE = 1000

# Proportion of full-length HERV
SIMULATE_RNA_VARIANTS_HERV_PROPORTION_FULL_LENGTH = 0.052 # She et al., Genome Biology 2022

# Enforce infinite sites assumption
SIMULATE_RNA_VARIANTS_INFINITE_SITES_ASSUMPTION = True


"""sim-reads"""
# Mean value of read length
SIMULATE_READS_READ_LENGTH_MEAN = 5000 # 9.759

# Standard deviation of read length
SIMULATE_READS_READ_LENGTH_STDEV = 500

# Mean value of base quality
SIMULATE_READS_BASE_QUALITY_MEAN = 90

# Standard deviation of base quality
SIMULATE_READS_BASE_QUALITY_STDEV = 5

# gzip
SIMULATE_READS_GZIP = True


"""identify"""
# Number of cores
NUM_CORES = 4


"""sim-meiosis"""
# Number of meitotic divisions for simulation of meiosis
NUM_MEITOTIC_DIVISIONS = 5

# Number of gametes to sample
NUM_SAMPLE_GAMETES = 100
