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
GENOME_GAPPED_REGIONS_PADDING = 1E5

# Padding for an excluded structural variant breakpoint.
# The pad is applied to upstream and downstream of the two breakpoints
# of a structural variant to be excluded.
EXCLUDE_SV_PADDING = 20

# Padding for an excluded single-nucleotide variant or
# small insertion / deletion (INDEL).
EXCLUDE_SNV_INDEL_PADDING = 1

# Maximum distance (bases) for merging two structural variants.
MAX_SV_CLUSTER_DISTANCE = 100

