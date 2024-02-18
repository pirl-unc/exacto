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


NUM_THREADS = 4


"""call-dna-vars"""
CALL_DNA_VARS_MIN_READS = 1
CALL_DNA_VARS_MIN_MAPPING_QUALITY = 0
CALL_DNA_VARS_MIN_INS_SIZE_PROPORTION = 0.5
CALL_DNA_VARS_MAX_INS_NORM_EDIT_DISTANCE = 0.5
CALL_DNA_VARS_MIN_DEL_SIZE_PROPORTION = 0.5
CALL_DNA_VARS_MAX_BND_DISTANCE = 1000
CALL_DNA_VARS_CLUSTERING_GRID_SIZE = 1000


"""call-rna-vars"""
CALL_RNA_VARS_MIN_READS = 3
CALL_RNA_VARS_MIN_MAPPING_QUALITY = 20
CALL_RNA_VARS_MIN_INS_SIZE_PROPORTION = 0.5
CALL_RNA_VARS_MAX_INS_NORM_EDIT_DISTANCE = 0.5
CALL_RNA_VARS_MIN_DEL_SIZE_PROPORTION = 0.5
CALL_RNA_VARS_MAX_BND_DISTANCE = 1000
CALL_RNA_VARS_CLUSTERING_GRID_SIZE = 1000


"""call-pep-vars"""
CALL_PEP_VARS_MIN_K = 9
CALL_PEP_VARS_MAX_K = 18
