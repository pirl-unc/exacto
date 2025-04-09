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

use crate::structs::vargraph_segment::VarGraphSegment;


#[derive(Clone,Debug)]
pub struct VarGraphPath {
    pub segments: Vec<VarGraphSegment>
}

impl Hash for VarGraphPath {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.segments.hash(state);
    }
}

impl PartialEq for VarGraphPath {
    fn eq(&self, other: &Self) -> bool {
        for i in 0..self.segments.len() {
            if self.segments[i] != other.segments[i] {
                return false;
            }
        }
        true
    }
}

impl VarGraphPath {
    pub fn new() -> Self {
        VarGraphPath {
            segments: Vec::new()
        }
    }

    pub fn display(&self) {
        println!("-------------------------------------------------------------------------");
        println!("VarGraphPath                                                             ");
        println!("-------------------------------------------------------------------------");
        println!("Node ID   | Strand   | CS tag                                           |");
        println!("--------- | -------- | ------------------------------------------------ |");
        for segment in self.segments.iter() {
            println!("{:9} | {:8} | {:48} |",
                     segment.node_id,
                     segment.strand.as_str(),
                     segment.cs_tag);
        }
        println!("-------------------------------------------------------------------------");
    }

    pub fn get_node_ids(&self) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        for segment in &self.segments {
            node_ids.push(segment.node_id);
        }
        node_ids
    }

    pub fn add_segment(&mut self, segment: &VarGraphSegment) {
        self.segments.push(segment.clone());
    }
}
