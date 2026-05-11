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

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraphSegment {
    pub node_id: usize,
    pub node_start: u32,
    pub node_end: u32,
    pub node_subsequence: Box<str>,
    pub enabled: bool,
    pub prev_node_id: Option<usize>,
    pub next_node_id: Option<usize>,
    pub entry_orientation: Option<VarGraphOrientations>,
    pub exit_orientation: Option<VarGraphOrientations>
}

impl PartialEq for VarGraphSegment {
    fn eq(&self, other: &Self) -> bool {
        self.node_id == other.node_id &&
            self.node_start == other.node_start &&
            self.node_end == other.node_end &&
            self.node_subsequence == other.node_subsequence &&
            self.enabled == other.enabled &&
            self.prev_node_id == other.prev_node_id &&
            self.next_node_id == other.next_node_id &&
            self.entry_orientation == other.entry_orientation &&
            self.exit_orientation == other.exit_orientation
    }
}

impl Hash for VarGraphSegment {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node_id.hash(state);
        self.node_start.hash(state);
        self.node_end.hash(state);
        self.node_subsequence.hash(state);
        self.enabled.hash(state);
        self.prev_node_id.hash(state);
        self.next_node_id.hash(state);
        self.entry_orientation.hash(state);
        self.exit_orientation.hash(state);
    }
}

impl VarGraphSegment {
    pub fn new(
        node_id: usize,
        node_start: u32,
        node_end: u32,
        node_subsequence: &str,
        prev_node_id: Option<usize>,
        next_node_id: Option<usize>,
        entry_orientation: Option<VarGraphOrientations>,
        exit_orientation: Option<VarGraphOrientations>
    ) -> Self {
        assert!(node_start <= node_end);
        VarGraphSegment {
            node_id: node_id,
            node_start: node_start,
            node_end: node_end,
            node_subsequence: node_subsequence.into(),
            enabled: true,
            prev_node_id: prev_node_id,
            next_node_id: next_node_id,
            entry_orientation: entry_orientation,
            exit_orientation: exit_orientation
        }
    }

    pub fn disable(&mut self) {
        self.enabled = false;
    }

    pub fn enable(&mut self) {
        self.enabled = true;
    }

    pub fn is_enabled(&self) -> bool {
        self.enabled
    }
}

impl Clone for VarGraphSegment {
    fn clone(&self) -> VarGraphSegment {
        let mut segment: VarGraphSegment = VarGraphSegment::new(
            self.node_id,
            self.node_start,
            self.node_end,
            &*self.node_subsequence,
            self.prev_node_id.clone(),
            self.next_node_id.clone(),
            self.entry_orientation.clone(),
            self.exit_orientation.clone()
        );
        segment.enabled = self.enabled;
        segment
    }
}
