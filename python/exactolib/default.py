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


"""call-dna-vars"""
CALL_DNA_VARS_MIN_READS = 3
CALL_DNA_VARS_MIN_MAPPING_QUALITY = 20
CALL_DNA_VARS_MIN_AVERAGE_BASE_QUALITY = 20
CALL_DNA_VARS_MIN_SIZE_PROPORTION = 0.5
CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE = 0.5
CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE = 1000
CALL_DNA_VARS_MAX_INTRACHROMOSOMAL_DISTANCE_TAU = 2000
CALL_DNA_VARS_MAX_INTERCHROMOSOMAL_DISTANCE = 1000
CALL_DNA_VARS_INFINITE_SITES_ASSUMPTION = 'yes'
CALL_DNA_VARS_NUM_THREADS = 4
CALL_DNA_VARS_GZIP = 'yes'

"""call-rna-vars"""
CALL_RNA_VARS_MIN_READS = 3
CALL_RNA_VARS_MIN_MAPPING_QUALITY = 20
CALL_RNA_VARS_MIN_AVERAGE_BASE_QUALITY = 20
CALL_RNA_VARS_NUM_THREADS = 4
CALL_RNA_VARS_REFERENCE_TRANSCRIPT_SCORING_METHOD = "weighted_net_overlap"

"""call-peptide-vars"""
CALL_PEPTIDE_VARS_K = 8
CALL_PEPTIDE_VARS_NUM_THREADS = 4
CALL_PEPTIDE_VARS_MIN_READS = 3
CALL_PEPTIDE_VARS_GZIP = 'yes'
CALL_PEPTIDE_VARS_DNA_VARIANT_PADDING = 100000

"""translate"""
TRANSLATE_STRATEGY = 'longest-orf'
TRANSLATE_NUM_THREADS = 4
TRANSLATE_GZIP = 'yes'
