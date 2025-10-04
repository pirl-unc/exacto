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


#[derive(Debug)]
pub struct VarGraphEdge {
    pub from: VarGraphPort,
    pub to: VarGraphPort,
    pub enabled: bool
}

impl PartialEq for VarGraphEdge {
    fn eq(&self, other: &Self) -> bool {
        self.from == other.from &&
            self.to == other.to &&
            self.enabled == other.enabled
    }
}

impl Hash for VarGraphEdge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.from.hash(state);
        self.to.hash(state);
        self.enabled.hash(state);
    }
}

impl VarGraphEdge {
    pub fn new(
        from: VarGraphPort,
        to: VarGraphPort
    ) -> Self {
        VarGraphEdge {
            from: from,
            to: to,
            enabled: true
        }
    }
}

impl Clone for VarGraphEdge {
    fn clone(&self) -> VarGraphEdge {
        let mut edge: VarGraphEdge = VarGraphEdge::new(
            self.from.clone(),
            self.to.clone()
        );
        edge.enabled = self.enabled;
        edge
    }
}
