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
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct AlignmentStructureEvent {
    prev_read_position: u32,
    next_read_position: u32,
    prev_graph_operation_type: GraphOperationType,
    next_graph_operation_type: GraphOperationType,
    kind: AlignmentStructureEventKind,
    context: Option<AlignmentStructureEventContext>,

    // For each chromosome ID (key), the vector (value) stores the 
    // reference positions in ascending order
    skipped_reference_bases: HashMap<u16, Vec<ReferenceBase>>
}

impl Hash for AlignmentStructureEvent {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kind.hash(state);
        self.prev_read_position.hash(state);
        self.next_read_position.hash(state);
        self.prev_graph_operation_type.hash(state);
        self.next_graph_operation_type.hash(state);
        self.context.hash(state);
    }
}

impl PartialEq for AlignmentStructureEvent {
    fn eq(&self, other: &Self) -> bool {
        self.kind == other.kind &&
            self.prev_read_position == other.prev_read_position &&
            self.next_read_position == other.next_read_position &&
            self.prev_graph_operation_type == other.prev_graph_operation_type &&
            self.next_graph_operation_type == other.next_graph_operation_type &&
            self.context == other.context
    }
}

impl AlignmentStructureEvent {
    pub fn new(
        kind: AlignmentStructureEventKind,
        prev_read_position: u32,
        next_read_position: u32,
        prev_graph_operation_type: GraphOperationType,
        next_graph_operation_type: GraphOperationType
    ) -> Self {
        assert!(
            prev_read_position <= next_read_position,
            "prev_read_position must be equal to or smaller than next_read_position."
        );
        Self {
            kind,
            prev_read_position,
            next_read_position,
            prev_graph_operation_type: prev_graph_operation_type,
            next_graph_operation_type: next_graph_operation_type,
            context: None,
            skipped_reference_bases: HashMap::new()
        }
    }

    pub fn add_reference_base(&mut self, reference_base: ReferenceBase) {
        // Get or create the chromosome vector
        let vec: &mut Vec<ReferenceBase> = self
            .skipped_reference_bases
            .entry(reference_base.reference_chromosome_id)
            .or_insert_with(Vec::new);

        // Enforce ascending position order
        let idx = vec
            .binary_search_by_key(&reference_base.reference_position, |b| b.reference_position)
            .unwrap_or_else(|idx| idx);
        vec.insert(idx, reference_base);
    }

    pub fn get_context(&self) -> Option<&AlignmentStructureEventContext> {
        self.context.as_ref()
    }
    
    pub fn get_kind(&self) -> &AlignmentStructureEventKind {
        &self.kind
    }

    pub fn get_prev_read_position(&self) -> u32 {
        self.prev_read_position
    }

    pub fn get_prev_graph_operation_type(&self) -> &GraphOperationType {
        &self.prev_graph_operation_type
    }
    
    pub fn get_next_read_position(&self) -> u32 {
        self.next_read_position
    }
    
    pub fn get_next_graph_operation_type(&self) -> &GraphOperationType {
        &self.next_graph_operation_type
    }

    pub fn get_skipped_reference_bases(&self) -> &HashMap<u16, Vec<ReferenceBase>> {
        &self.skipped_reference_bases
    }

    pub fn get_skipped_reference_bases_clusters(&self) -> Vec<Vec<&ReferenceBase>> {
        let mut clusters: Vec<Vec<&ReferenceBase>> = Vec::new();
        for (chromosome_id, bases) in self.skipped_reference_bases.iter() {
            // Cluster reference bases by proximity
            let mut uf: UnionFind = UnionFind::new();
            for i in 0..bases.len() {
                uf.union(i, i);
                if i > 0 {
                    if bases[i].reference_position.abs_diff(bases[i - 1].reference_position) == 1 {
                        uf.union(i - 1, i);
                    }
                }
            }

            // Get the clusters
            for cluster in uf.get_clusters() {
                // Sort the cluster
                let mut indices: Vec<usize> = cluster.into_iter().collect();
                indices.sort();

                // Get the cluster bases
                let selected: Vec<&ReferenceBase> = indices.iter()
                    .map(|&i| bases.get(i).unwrap())
                    .collect();

                clusters.push(selected);
            }
        }

        clusters
    }

    pub fn set_context(&mut self, context: AlignmentStructureEventContext) {
        self.context = Some(context);
    }

    pub fn set_kind(&mut self, kind: AlignmentStructureEventKind) {
        self.kind = kind;
    }
}

impl Clone for AlignmentStructureEvent {
    fn clone(&self) -> Self {
        AlignmentStructureEvent {
            kind: self.kind.clone(),
            prev_read_position: self.prev_read_position,
            next_read_position: self.next_read_position,
            prev_graph_operation_type: self.prev_graph_operation_type.clone(),
            next_graph_operation_type: self.next_graph_operation_type.clone(),
            context: self.context.clone(),
            skipped_reference_bases: self.skipped_reference_bases.clone()
        }
    }
}
