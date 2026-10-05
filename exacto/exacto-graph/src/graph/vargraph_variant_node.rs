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


use exacto_caller::prelude::*;
use exacto_core::prelude::*;

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraphVariantNode {
    pub variant_id: usize,
    pub graph_operation_view: GraphOperationView
}

impl VarGraphVariantNode {
    pub fn new(
        variant_id: usize,
        graph_operation_view: GraphOperationView
    ) -> Self {
        VarGraphVariantNode {
            variant_id: variant_id,
            graph_operation_view: graph_operation_view
        }
    }

    pub fn get_chromosome_1(&self) -> &str {
        &self.graph_operation_view.get_chromosome_1()
    }
    
    pub fn get_chromosome_2(&self) -> &str {
        &self.graph_operation_view.get_chromosome_2()
    }
    
    pub fn get_position_1(&self) -> u32 {
        self.graph_operation_view.get_position_1()
    }
    
    pub fn get_position_2(&self) -> u32 {
        self.graph_operation_view.get_position_2()
    }
    
    pub fn get_operation_1(&self) -> &GraphOperationType {
        self.graph_operation_view.get_operation_type_1()
    }

    pub fn get_operation_2(&self) -> &GraphOperationType {
        self.graph_operation_view.get_operation_type_2()
    }
    
    pub fn get_strand_1(&self) -> &Strand {
        self.graph_operation_view.get_strand_1()
    }
    
    pub fn get_strand_2(&self) -> &Strand {
        self.graph_operation_view.get_strand_2()
    }
    
    /// The node's bases, read from its position_1 side to its position_2 side.
    ///
    /// The variant tables spell a sequence genome-forward, in the order a forward read passes the
    /// junction (exacto-caller writes `get_standardized_sequence`), so the strand does not change
    /// it. A forward read passes a D-to-U junction from position_1 to position_2, but a U-to-D
    /// junction (a tandem duplication, or a translocation joined that way) from position_2 to
    /// position_1. A path walks such a node from its position_2 side and reverse-complements it,
    /// so the node holds the reverse complement and the path writes the bases as stored.
    pub fn get_sequence(&self) -> Box<str> {
        let sequence: String = self.graph_operation_view.get_sequence().to_uppercase();
        if *self.get_operation_1() == GraphOperationType::Upstream &&
            *self.get_operation_2() == GraphOperationType::Downstream {
            reverse_complement(&sequence)
        } else {
            sequence.into_boxed_str()
        }
    }

    pub fn get_sequence_length(&self) -> u32 {
        self.graph_operation_view.get_sequence_length() as u32
    }
    
    pub fn get_variant_id(&self) -> usize {
        self.variant_id
    }
}

impl Clone for VarGraphVariantNode {
    fn clone(&self) -> VarGraphVariantNode {
        VarGraphVariantNode::new(
            self.variant_id,
            self.graph_operation_view.clone()
        )
    }
}
