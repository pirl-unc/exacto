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

# Probability of simulating a genic variant
SIMULATE_GENIC_VARIANT_PROBABILITY = 0.2

# Number of single-nucleotide variants to simulate
SIMULATE_NUM_SNV = 100

# Number of small insertions to simulate
SIMULATE_NUM_INSERTION = 100

# Number of small deletions to simulate
SIMULATE_NUM_DELETION = 100

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
    'genotype': 'unknown',                                                     # genotype
    'genotype_quality': 'unknown',                                             # genotype quality
    'total_coverage': 'unknown',                                               # total coverage
    'reference_reads_count': 'unknown',                                        # reference reads count
    'variant_reads_count': 'unknown',                                          # variant reads count
    'variant_allele_fraction': 'unknown',                                      # variant allele fraction
    'phred_scale_genotype_likelihoods': 'unknown'
}
