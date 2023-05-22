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


"""convert-variant-to-tsv"""


"""merge-variant-calls"""
# Maximum neighbor distance (bases).
MERGE_MAX_NEIGHBOR_DISTANCE = 1

# Number of processes
MERGE_NUM_PROCESSES = 4


"""filter-variants"""
# Padding for an excluded region.
# The pad is applied to upstream and downstream of a gapped genomic region.
FILTER_VARIANTS_EXCLUDED_REGION_PADDING = 100000

# Padding for an excluded variant.
# The pad is applied to upstream and downstream of the two breakpoints
# of a variant to be excluded.
FILTER_VARIANTS_EXCLUDED_VARIANT_PADDING = 1000

# Number of processes
FILTER_VARIANTS_NUM_PROCESSES = 4


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


"""sim-variants"""
# Probability of simulating a genic variant
SIMULATE_VARIANTS_GENIC_PROBABILITY = 0.2

# Number of single-nucleotide variants to simulate
SIMULATE_VARIANTS_NUM_SNV = 300

# Number of insertions to simulate
SIMULATE_VARIANTS_NUM_INSERTION = 50

# Insertion size mean
SIMULATE_VARIANTS_INSERTION_SIZE_MEAN = 100

# Insertion size standard deviation
SIMULATE_VARIANTS_INSERTION_SIZE_STDEV = 50

# Number of deletions to simulate
SIMULATE_VARIANTS_NUM_DELETION = 50

# Deletion size mean
SIMULATE_VARIANTS_DELETION_MEAN = 100

# Deletion size standard deviation
SIMULATE_VARIANTS_DELETION_STDEV = 50

# Enforce infinite sites model
SIMULATE_VARIANTS_ENFORCE_INFINITE_SITES_MODEL = True


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


"""call-variants"""
# Number of processes
CALL_VARIANTS_NUM_PROCESSES = 4


"""sim-meiosis"""
# Number of meitotic divisions for simulation of meiosis
NUM_MEITOTIC_DIVISIONS = 5

# Number of gametes to sample
NUM_SAMPLE_GAMETES = 100
