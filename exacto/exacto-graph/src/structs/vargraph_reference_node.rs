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


use exacto_util::prelude::*;
use std::any::Any;

use crate::common::constants::VarGraphNodeTypes;
use crate::traits::vargraph_node::VarGraphNode;


#[derive(Debug)]
pub struct VarGraphReferenceNode {
    pub chromosome: Box<str>,
    pub start: usize,
    pub end: usize,
    pub sequence: Box<str>
}

impl VarGraphReferenceNode {
    pub fn new(chromosome: &str, start: usize, end: usize, sequence: &str) -> Self {
        VarGraphReferenceNode {
            chromosome: chromosome.into(),
            start: start,
            end: end,
            sequence: sequence.into()
        }
    }
}

impl VarGraphNode for VarGraphReferenceNode {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn clone_box(&self) -> Box<dyn VarGraphNode> {
        Box::new((*self).clone())
    }

    fn get_chromosome_1(&self) -> Box<str> {
        self.chromosome.clone()
    }

    fn get_chromosome_2(&self) -> Box<str> {
        self.chromosome.clone()
    }

    fn get_position_1(&self) -> usize {
        self.start
    }

    fn get_position_2(&self) -> usize {
        self.end
    }

    fn get_reverse_complement_sequence(&self) -> Box<str> {
        reverse_complement(&*self.sequence)
    }

    fn get_sequence(&self) -> Box<str> {
        self.sequence.clone()
    }

    fn get_sequence_length(&self) -> isize {
        (self.end - self.start + 1) as isize
    }

    fn get_type(&self) -> VarGraphNodeTypes {
        VarGraphNodeTypes::Reference
    }
}

impl Clone for VarGraphReferenceNode {
    fn clone(&self) -> VarGraphReferenceNode {
        let reference_node: VarGraphReferenceNode = VarGraphReferenceNode::new(
            &*self.chromosome,
            self.start,
            self.end,
            &*self.sequence
        );
        reference_node
    }
}
