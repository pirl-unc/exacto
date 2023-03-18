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


"""merge"""
# Maximum neighbor distance (bases).
MAX_NEIGHBOR_DISTANCE = 1


"""filter"""
# Padding for an excluded region.
# The pad is applied to upstream and downstream of a gapped genomic region.
EXCLUDED_REGION_PADDING = 100000

# Padding for an excluded variant.
# The pad is applied to upstream and downstream of the two breakpoints
# of a variant to be excluded.
EXCLUDED_VARIANT_PADDING = 100

# Enforce variant type matching.
ENFORCE_VARIANT_TYPE_MATCHING = True

# Number of processes
NUM_PROCESSES_REFINE = 4


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
