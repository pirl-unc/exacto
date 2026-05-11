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


use exacto_core::prelude::*;
use std::collections::VecDeque;
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Clone,Debug)]
pub struct VarGraphPath {
    pub segments: VecDeque<VarGraphSegment>
}

impl Hash for VarGraphPath {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.segments.hash(state);
    }
}

impl PartialEq for VarGraphPath {
    fn eq(&self, other: &Self) -> bool {
        if self.segments.len() != other.segments.len() {
            return false;
        }
        for i in 0..self.segments.len() {
            if self.segments[i] != other.segments[i] {
                return false;
            }
        }
        true
    }
}

impl Eq for VarGraphPath {}

/// API methods
impl VarGraphPath {
    pub fn new() -> Self {
        VarGraphPath {
            segments: VecDeque::new()
        }
    }

    pub fn push_back(&mut self, segment: VarGraphSegment) {
        self.segments.push_back(segment);
    }
    
    pub fn push_front(&mut self, segment: VarGraphSegment) {
        self.segments.push_front(segment);
    }

    pub fn get_first_segment(&self) -> &VarGraphSegment {
        self.segments.front().unwrap()
    }

    pub fn get_last_segment(&self) -> &VarGraphSegment {
        self.segments.back().unwrap()
    }

    pub fn get_node_ids(&self) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        for segment in &self.segments {
            node_ids.push(segment.node_id);
        }
        node_ids
    }

    pub fn get_sequence(&self) -> Box<str> {
        let mut sequence: String = String::new();
        for (i, segment) in self.segments.iter().enumerate() {
            if segment.is_enabled() == false {
                continue;
            }
            let should_rc = if let Some(entry) = segment.entry_orientation.as_ref() {
                *entry != VarGraphOrientations::Upstream
            } else if i == 0 {
                // Genome paths: first segment has no entry_orientation,
                // fall back to exit_orientation
                segment.exit_orientation.as_ref()
                    .map(|exit| *exit == VarGraphOrientations::Upstream)
                    .unwrap_or(false)
            } else {
                false
            };
            if should_rc {
                sequence.push_str(&*reverse_complement(&*segment.node_subsequence));
            } else {
                sequence.push_str(&*segment.node_subsequence);
            }
        }

        sequence.into_boxed_str()
    }
}
