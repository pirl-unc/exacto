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


#[derive(Debug, Clone)]
pub struct IdentifyConsensusSequencesOptions {
    pub poa_match_score: i32,
    pub poa_mismatch_score: i32,
    pub poa_gap_open_score: i32,
    pub poa_gap_extend_score: i32,

    /// A cluster of more reads than this is aligned on this many of its reads, drawn at random.
    /// 0 aligns all the reads of every cluster.
    pub max_reads_per_cluster: usize,

    /// Seed for the subsamples. The seed of a cluster is this seed plus its cluster ID.
    pub seed: usize,

    /// Length of the k-mers by which the reads of a cluster are turned to the orientation of
    /// most of them before alignment. At most 32. 0 aligns the reads as they are in the FASTQ.
    pub orientation_kmer_size: usize
}

/// IdentifyConsensusSequencesOptions
impl IdentifyConsensusSequencesOptions {
    pub const DEFAULT: Self = Self {
        poa_match_score: 0,
        poa_mismatch_score: 4,
        poa_gap_open_score: 6,
        poa_gap_extend_score: 2,
        max_reads_per_cluster: 1_000,
        seed: 42,
        orientation_kmer_size: 15
    };
}
impl Default for IdentifyConsensusSequencesOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}