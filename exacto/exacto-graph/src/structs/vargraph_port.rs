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

use crate::common::constants::{VarGraphOrientations, VarGraphStrands};


#[derive(Debug)]
pub struct VarGraphPort {
    pub node_id: usize,
    pub strand: VarGraphStrands,
    pub orientation: VarGraphOrientations
}

impl PartialEq for VarGraphPort {
    fn eq(&self, other: &Self) -> bool {
        self.node_id == other.node_id &&
            self.strand == other.strand &&
            self.orientation == other.orientation
    }
}

impl Hash for VarGraphPort {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node_id.hash(state);
        self.strand.as_str().hash(state);
        self.orientation.as_str().hash(state);
    }
}

impl VarGraphPort {
    pub fn new(
        node_id: usize,
        strand: VarGraphStrands,
        orientation: VarGraphOrientations
    ) -> Self {
        VarGraphPort {
            node_id: node_id,
            strand: strand,
            orientation: orientation
        }
    }
}

impl Clone for VarGraphPort {
    fn clone(&self) -> VarGraphPort {
        let port: VarGraphPort = VarGraphPort::new(
            self.node_id,
            self.strand.clone(),
            self.orientation.clone()
        );
        port
    }
}
