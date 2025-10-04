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
use std::any::Any;

use crate::common::constants::*;
use crate::traits::vargraph_node::VarGraphNode;


#[derive(Debug)]
pub struct VarGraphVariantNode {
    pub chromosome_1: Box<str>,
    pub position_1: usize,
    pub orientation_1: VarGraphOrientations,
    pub strand_1: VarGraphStrands,
    pub chromosome_2: Box<str>,
    pub position_2: usize,
    pub orientation_2: VarGraphOrientations,
    pub strand_2: VarGraphStrands,
    pub sequence: Box<str>
}

impl VarGraphVariantNode {
    pub fn new(
        chromosome_1: &str,
        position_1: usize,
        orientation_1: VarGraphOrientations,
        strand_1: VarGraphStrands,
        chromosome_2: &str,
        position_2: usize,
        orientation_2: VarGraphOrientations,
        strand_2: VarGraphStrands,
        sequence: &str
    ) -> Self {
        VarGraphVariantNode {
            chromosome_1: chromosome_1.into(),
            position_1: position_1,
            orientation_1: orientation_1,
            strand_1: strand_1,
            chromosome_2: chromosome_2.into(),
            position_2: position_2,
            orientation_2: orientation_2,
            strand_2: strand_2,
            sequence: sequence.into()
        }
    }
}

impl VarGraphNode for VarGraphVariantNode {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn clone_box(&self) -> Box<dyn VarGraphNode> {
        Box::new((*self).clone())
    }

    fn get_chromosome_1(&self) -> Box<str> {
        self.chromosome_1.clone()
    }

    fn get_chromosome_2(&self) -> Box<str> {
        self.chromosome_2.clone()
    }

    fn get_position_1(&self) -> usize {
        self.position_1
    }

    fn get_position_2(&self) -> usize {
        self.position_2
    }

    fn get_reverse_complement_sequence(&self) -> Box<str> {
        if self.sequence != "".into() {
            reverse_complement(&*self.sequence)
        } else {
            "".to_string().into_boxed_str()
        }
    }

    fn get_sequence(&self) -> Box<str> {
        self.sequence.clone()
    }

    fn get_sequence_length(&self) -> isize {
        self.sequence.len() as isize
    }

    fn get_type(&self) -> VarGraphNodeTypes {
        VarGraphNodeTypes::Variant
    }
}

impl Clone for VarGraphVariantNode {
    fn clone(&self) -> VarGraphVariantNode {
        let node: VarGraphVariantNode = VarGraphVariantNode::new(
            &*self.chromosome_1,
            self.position_1,
            self.orientation_1.clone(),
            self.strand_1.clone(),
            &*self.chromosome_2,
            self.position_2,
            self.orientation_2.clone(),
            self.strand_2.clone(),
            &*self.sequence
        );
        node
    }
}
