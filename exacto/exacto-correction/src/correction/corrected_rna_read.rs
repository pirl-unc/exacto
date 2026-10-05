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


use serde::{Serialize, Deserialize};
use std::hash::{Hash, Hasher};


#[derive(Debug, Default, Serialize, Deserialize, Clone)]
pub struct CorrectedRNARead {
    pub cluster_id: usize,
    pub read_name: Box<str>,
    pub sequence: Box<str>,
    pub base_quality_scores: Vec<u8>
}

impl PartialEq for CorrectedRNARead {
    fn eq(&self, other: &Self) -> bool {
        self.cluster_id == other.cluster_id 
            && self.read_name == other.read_name
            && self.sequence == other.sequence
            && self.base_quality_scores == other.base_quality_scores
    }
}

impl Eq for CorrectedRNARead {}

impl Hash for CorrectedRNARead {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cluster_id.hash(state);
        self.read_name.hash(state);
        self.sequence.hash(state);
        self.base_quality_scores.hash(state);
    }
}

impl CorrectedRNARead {
    pub fn new(
        cluster_id: usize,
        read_name: Box<str>,
        sequence: Box<str>,
        base_quality_scores: Vec<u8>
    ) -> Self {
        Self {
            cluster_id: cluster_id,
            read_name: read_name,
            sequence: sequence,
            base_quality_scores: base_quality_scores
        }
    }
}
