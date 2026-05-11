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


use crate::prelude::*;


#[derive(Debug, Clone)]
pub enum VarGraphNode {
    Reference(VarGraphReferenceNode),
    Variant(VarGraphVariantNode)
}

impl VarGraphNode {
    pub fn as_reference(&self) -> &VarGraphReferenceNode {
        match self {
            Self::Reference(node) => node,
            _ => panic!("Not a reference node.")
        }
    }

    pub fn as_reference_mut(&mut self) -> &mut VarGraphReferenceNode {
        match self {
            Self::Reference(node) => node,
            _ => panic!("Not a reference node."),
        }
    }
    
    pub fn as_variant(&self) -> &VarGraphVariantNode {
        match self {
            Self::Variant(node) => node,
            _ => panic!("Not a variant node.")
        }
    }

    pub fn as_variant_mut(&mut self) -> &mut VarGraphVariantNode {
        match self {
            Self::Variant(node) => node,
            _ => panic!("Not a variant node.")
        }
    }
    
    pub fn get_chromosome_1(&self) -> &str {
        match self {
            Self::Reference(node) => {
                &node.get_chromosome()
            },
            Self::Variant(node) => {
                &node.graph_operation_view.get_chromosome_1()
            }
        }
    }

    pub fn get_chromosome_2(&self) -> &str {
        match self {
            Self::Reference(node) => {
                &node.get_chromosome()
            },
            Self::Variant(node) => {
                &node.graph_operation_view.get_chromosome_2()
            }
        }
    }

    pub fn get_position_1(&self) -> u32 {
        match self {
            Self::Reference(node) => {
                node.get_start()
            },
            Self::Variant(node) => {
                node.graph_operation_view.get_position_1()
            }
        }
    }

    pub fn get_position_2(&self) -> u32 {
        match self {
            Self::Reference(node) => {
                node.get_end()
            },
            Self::Variant(node) => {
                node.graph_operation_view.get_position_2()
            }
        }
    }

    pub fn get_sequence(&self) -> Box<str> {
        match self {
            Self::Reference(node) => {
                node.get_sequence().into()
            },
            Self::Variant(node) => {
                node.graph_operation_view.get_standardized_sequence()
            }
        }
    }

    pub fn get_sequence_length(&self) -> u32 {
        match self {
            Self::Reference(node) => {
                node.get_sequence_length()
            },
            Self::Variant(node) => {
                node.graph_operation_view.get_standardized_sequence().len() as u32
            }
        }
    }

    pub fn get_type(&self) -> VarGraphNodeTypes {
        match self {
            Self::Reference(_) => VarGraphNodeTypes::Reference,
            Self::Variant(_) => VarGraphNodeTypes::Variant
        }
    }
}
