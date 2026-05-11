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
    
    pub fn get_sequence(&self) -> Box<str> {
        self.graph_operation_view.get_standardized_sequence()
    }
    
    pub fn get_sequence_length(&self) -> u32 {
        self.graph_operation_view.get_standardized_sequence().len() as u32
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
