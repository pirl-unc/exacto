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


use std::collections::{HashSet, VecDeque};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraphTraversalState {
    pub visited: HashSet<usize>,
    pub path: VecDeque<usize>,
    pub prev_edge: Option<VarGraphEdge>
}

impl Hash for VarGraphTraversalState {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.path.hash(state);
    }
}

impl PartialEq for VarGraphTraversalState {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl Eq for VarGraphTraversalState {}

impl VarGraphTraversalState {
    pub fn new(
        visited: HashSet<usize>,
        path: VecDeque<usize>,
        prev_edge: Option<VarGraphEdge>
    ) -> Self {
        Self {
            visited: visited,
            path: path,
            prev_edge: prev_edge
        }
    }
}

impl Clone for VarGraphTraversalState {
    fn clone(&self) -> Self {
        VarGraphTraversalState {
            visited: self.visited.clone(),
            path: self.path.clone(),
            prev_edge: self.prev_edge.clone()
        }
    }
}
