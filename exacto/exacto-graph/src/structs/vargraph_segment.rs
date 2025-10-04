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


use std::hash::{Hash, Hasher};

use crate::structs::vargraph_port::VarGraphPort;


#[derive(Clone,Debug)]
pub struct VarGraphSegment {
    pub node_id: usize,
    pub node_start: usize,
    pub node_end: usize,
    pub node_subsequence: Box<str>,             // in the strand of entry_port.strand
    pub entry_port: Option<VarGraphPort>,
    pub exit_port: Option<VarGraphPort>
}

impl PartialEq for VarGraphSegment {
    fn eq(&self, other: &Self) -> bool {
        self.node_id == other.node_id &&
            self.node_start == other.node_start &&
            self.node_end == other.node_end &&
            self.node_subsequence == other.node_subsequence &&
            self.entry_port == other.entry_port &&
            self.exit_port == other.exit_port
    }
}

impl Hash for VarGraphSegment {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node_id.hash(state);
        self.node_start.hash(state);
        self.node_end.hash(state);
        self.node_subsequence.hash(state);
        self.entry_port.hash(state);
    }
}

impl VarGraphSegment {
    pub fn new(
        node_id: usize,
        node_start: usize,
        node_end: usize,
        node_subsequence: &str,
        entry_port: Option<VarGraphPort>,
        exit_port: Option<VarGraphPort>
    ) -> Self {
        assert!(node_start <= node_end);
        VarGraphSegment {
            node_id: node_id,
            node_start: node_start,
            node_end: node_end,
            node_subsequence: node_subsequence.into(),
            entry_port: entry_port,
            exit_port: exit_port
        }
    }
}
