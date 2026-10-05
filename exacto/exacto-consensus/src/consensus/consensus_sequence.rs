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
use std::collections::HashSet;
use std::hash::{Hash, Hasher};


#[derive(Debug,Eq,Serialize,Deserialize)]
pub struct ConsensusSequence {
    /// Globally-unique id.
    pub cluster_id: usize,

    /// Consensus sequence.
    pub consensus_sequence: Box<str>,
    
    /// Reads in this cluster.
    pub read_names: HashSet<Box<str>>
}

impl PartialEq for ConsensusSequence {
    fn eq(&self, other: &Self) -> bool {
        if self.cluster_id == other.cluster_id {
            true
        } else {
            false
        }
    }
}

impl Hash for ConsensusSequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cluster_id.hash(state);
    }
}

impl ConsensusSequence {
    pub fn new(
        cluster_id: usize,
        consensus_sequence: Box<str>,
        read_names: HashSet<Box<str>>
    ) -> Self {
        Self {
            cluster_id: cluster_id,
            consensus_sequence: consensus_sequence,
            read_names: read_names
        }
    }

    pub fn get_cluster_id(&self) -> usize {
        self.cluster_id
    }

    pub fn get_consensus_sequence(&self) -> &str {
        &*self.consensus_sequence
    }

    pub fn get_read_names(&self) -> &HashSet<Box<str>> {
        &self.read_names
    }
}

impl Clone for ConsensusSequence {
    fn clone(&self) -> Self {
        ConsensusSequence {
            cluster_id: self.cluster_id,
            consensus_sequence: self.consensus_sequence.clone(),
            read_names: self.read_names.clone()
        }
    }
}
