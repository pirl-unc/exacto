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


use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use exacto_core::log_info;
use std::any::Any;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::ops::Bound::Included;
use std::str::FromStr;
use std::sync::Arc;

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraph {
    multidigraph: MultiDiGraph<VarGraphNode, VarGraphEdge>,

    /// Reference nodes index.
    ///
    /// HashMap
    ///     Key     :   Chromosome.
    ///     Value   :   BTreeMap (key = (start, end), value = node ID).
    reference_nodes_index: HashMap<Box<str>, BTreeMap<(u32, u32), usize>>,

    /// Variants
    ///
    /// HashMap
    ///     Key     :   Variant ID.
    ///     Value   :   GraphOperationView.
    operation_by_variant_id: HashMap<usize, GraphOperationView>,

    /// Variant IDs index.
    ///
    /// HashMap
    ///     Left    :   Variant ID.
    ///     Right   :   Node IDs.
    node_ids_by_variant_id: HashMap<usize, HashSet<usize>>,

    /// Variant operations index.
    ///
    /// HashMap:
    ///     Key     :   GraphOperationView.
    ///     Value   :   Node IDs.
    node_ids_by_operation: HashMap<GraphOperationView, HashSet<usize>>
}

/// API methods
impl VarGraph {
    /// Constructors
    pub fn new() -> Self{
        VarGraph {
            multidigraph: MultiDiGraph::new(),
            reference_nodes_index: HashMap::new(),
            operation_by_variant_id: HashMap::new(),
            node_ids_by_variant_id: HashMap::new(),
            node_ids_by_operation: HashMap::new()
        }
    }

    pub fn from_fasta_file(fasta_file: &str) -> Self {
        let mut vargraph = VarGraph {
            multidigraph: MultiDiGraph::new(),
            reference_nodes_index: HashMap::new(),
            operation_by_variant_id: HashMap::new(),
            node_ids_by_variant_id: HashMap::new(),
            node_ids_by_operation: HashMap::new()
        };

        let sequence_ids: Vec<(Box<str>, u32)> = get_fasta_sequence_ids(fasta_file);

        for (sequence_id, sequence_length) in &sequence_ids {
            let reference_sequence: Box<str> = get_fasta_sequence(&*sequence_id, 1, *sequence_length, fasta_file);
            let node: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*sequence_id,
                1,
                *sequence_length,
                &*reference_sequence
            );
            vargraph.add_reference_node(node);
        }

        vargraph
    }

    pub fn from_reference_nodes(reference_nodes_list: Vec<Vec<VarGraphReferenceNode>>) -> Self {
        let mut vargraph = VarGraph {
            multidigraph: MultiDiGraph::new(),
            reference_nodes_index: HashMap::new(),
            operation_by_variant_id: HashMap::new(),
            node_ids_by_variant_id: HashMap::new(),
            node_ids_by_operation: HashMap::new()
        };

        for reference_nodes in reference_nodes_list {
            let mut prev_reference_node_id: Option<usize> = None;
            for node in reference_nodes {
                let reference_node_id: usize = vargraph.add_reference_node(node);
                if prev_reference_node_id.is_some() {
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Downstream,
                        VarGraphOrientations::Upstream
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Upstream,
                        VarGraphOrientations::Downstream
                    );
                    vargraph.multidigraph.add_edge(prev_reference_node_id.unwrap(), reference_node_id, edge_1);
                    vargraph.multidigraph.add_edge(reference_node_id, prev_reference_node_id.unwrap(), edge_2);
                }
                prev_reference_node_id = Some(reference_node_id);
            }
        }

        vargraph
    }

    /// Adds a reference node.
    ///
    /// Parameters:
    /// - node          :   VarGraphReferenceNode object.
    ///
    /// Returns:
    /// - Reference node ID.
    pub fn add_reference_node(&mut self, node: VarGraphReferenceNode) -> usize {
        let chromosome: Box<str> = node.get_chromosome().into();
        let start: u32 = node.get_start();
        let end: u32 = node.get_end();

        // Step 1. Error checking
        assert!(start <= end, "start must be equal to or smaller than end");
        assert!(node.get_sequence().len() as u32 == end - start + 1, "sequence length must match the length from start to end");

        // Step 2. Add reference node
        let new_node_id: usize = self.multidigraph.add_node(VarGraphNode::Reference(node));

        // Step 3. Add the reference node ID to index
        self.reference_nodes_index
            .entry(chromosome.into())
            .or_insert_with(BTreeMap::new)
            .insert((start, end), new_node_id);

        new_node_id
    }

    /// Add a variant node.
    ///
    /// Parameters:
    /// - node      :   VarGraphVariantNode object.
    ///
    /// Returns:
    /// - Variant node ID.
    pub fn add_variant_node(&mut self, node: VarGraphVariantNode) -> Option<usize> {
        // Step 1. Check if the variant node already exists
        if let Some(existing_id) = self.query_variant_node(&node.graph_operation_view) {
            return Some(existing_id);
        }

        // Step 2. Get the node attributes
        let variant_id: usize = node.variant_id;
        let chromosome_1: Box<str> = node.graph_operation_view.get_chromosome_1().into();
        let position_1: u32 = node.graph_operation_view.get_position_1();
        let operation_type_1: GraphOperationType = node.graph_operation_view.get_operation_type_1().clone();
        let chromosome_2: Box<str> = node.graph_operation_view.get_chromosome_2().into();
        let position_2: u32 = node.graph_operation_view.get_position_2();
        let operation_type_2: GraphOperationType = node.graph_operation_view.get_operation_type_2().clone();
        let gov: GraphOperationView = node.graph_operation_view.clone();

        // Step 3. Add the variant node
        if (operation_type_1 == GraphOperationType::Downstream || operation_type_1 == GraphOperationType::Upstream) &&
            (operation_type_2 == GraphOperationType::Downstream || operation_type_2 == GraphOperationType::Upstream) {
            let new_node_id: usize = self.multidigraph.add_node(VarGraphNode::Variant(node));

            // Split the reference node at position_1 and position_2
            let mut split_reference_at_position = |chromosome: &str, position: u32| -> usize {
                let ref_node_id: usize = self.query_reference_node(chromosome, position)
                    .unwrap_or_else(|| {
                        panic!(
                            "Expected to have a reference node containing the position {}:{}",
                            chromosome, position
                        )
                    });
                self.split_reference_node(ref_node_id, position)
            };
            let pos_1_ref_node_id: usize = split_reference_at_position(&*chromosome_1, position_1);
            let pos_2_ref_node_id: usize = split_reference_at_position(&*chromosome_2, position_2);

            // Disable outgoing reference nodes from positions 1 and 2 unless the variant creates a cyclic path
            if !((chromosome_1 == chromosome_2) && (operation_type_1 == GraphOperationType::Upstream && operation_type_2 == GraphOperationType::Downstream)) {
                for outgoing_node_id in self.get_outgoing_node_ids(pos_1_ref_node_id).clone() {
                    let orientation_from: VarGraphOrientations = self
                        .get_edge_data(pos_1_ref_node_id, outgoing_node_id)
                        .orientation_from
                        .clone();
                    if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Reference {
                        if operation_type_1 == GraphOperationType::Downstream &&
                            orientation_from == VarGraphOrientations::Downstream {
                            self.disable_edge(pos_1_ref_node_id, outgoing_node_id);
                        }
                        if operation_type_1 == GraphOperationType::Upstream &&
                            orientation_from == VarGraphOrientations::Upstream {
                            self.disable_edge(pos_1_ref_node_id, outgoing_node_id);
                        }
                    }
                }
                for outgoing_node_id in self.get_outgoing_node_ids(pos_2_ref_node_id).clone() {
                    let orientation_from: VarGraphOrientations = self
                        .get_edge_data(pos_2_ref_node_id, outgoing_node_id)
                        .orientation_from
                        .clone();
                    if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Reference {
                        if operation_type_2 == GraphOperationType::Downstream &&
                            orientation_from == VarGraphOrientations::Downstream {
                            self.disable_edge(pos_2_ref_node_id, outgoing_node_id);
                        }
                        if operation_type_2 == GraphOperationType::Upstream &&
                            orientation_from == VarGraphOrientations::Upstream {
                            self.disable_edge(pos_2_ref_node_id, outgoing_node_id);
                        }
                    }
                }
            }

            // Rewrite position_1
            if operation_type_1 == GraphOperationType::Downstream {
                let edge_1: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Downstream,
                    VarGraphOrientations::Upstream
                );
                let edge_2: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Upstream,
                    VarGraphOrientations::Downstream
                );
                self.multidigraph.add_edge(pos_1_ref_node_id, new_node_id, edge_1);
                self.multidigraph.add_edge(new_node_id, pos_1_ref_node_id, edge_2);
            } else {
                let edge_1: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Upstream,
                    VarGraphOrientations::Upstream
                );
                let edge_2: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Upstream,
                    VarGraphOrientations::Upstream
                );
                self.multidigraph.add_edge(pos_1_ref_node_id, new_node_id, edge_1);
                self.multidigraph.add_edge(new_node_id, pos_1_ref_node_id, edge_2);
            }

            // Rewire position_2
            if operation_type_2 == GraphOperationType::Downstream {
                let edge_1: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Downstream,
                    VarGraphOrientations::Downstream
                );
                let edge_2: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Downstream,
                    VarGraphOrientations::Downstream
                );
                self.multidigraph.add_edge(new_node_id, pos_2_ref_node_id, edge_1);
                self.multidigraph.add_edge(pos_2_ref_node_id, new_node_id, edge_2);
            } else {
                let edge_1: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Downstream,
                    VarGraphOrientations::Upstream
                );
                let edge_2: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Upstream,
                    VarGraphOrientations::Downstream
                );
                self.multidigraph.add_edge(new_node_id, pos_2_ref_node_id, edge_1);
                self.multidigraph.add_edge(pos_2_ref_node_id, new_node_id, edge_2);
            }

            // Add the variant node ID to the indices
            self.operation_by_variant_id.insert(variant_id, gov.clone());
            self.node_ids_by_variant_id.entry(variant_id).or_insert_with(HashSet::new).insert(new_node_id);
            self.node_ids_by_operation.entry(gov.clone()).or_insert_with(HashSet::new).insert(new_node_id);

            return Some(new_node_id);
        }
        if operation_type_1 == GraphOperationType::Include && operation_type_2 == GraphOperationType::Include {
            assert!(chromosome_1 == chromosome_2);
            assert!(position_1 <= position_2);

            self.annotate_reference_nodes(
                &*chromosome_1,
                position_1,
                position_2,
                VarGraphReferenceNodeAnnotationKeys::Type.as_str(),
                Arc::new(VarGraphReferenceNodeTypes::Exonic.as_str().to_string().into_boxed_str()) as Arc<dyn Any + Send + Sync>
            );

            // Add the variant node ID to the indices
            self.operation_by_variant_id.insert(variant_id, gov.clone());
            let ref_node_ids: Vec<usize> = self.reference_nodes_index
                .get(&*chromosome_1)
                .unwrap()
                .range((Included((position_1, 0u32)), Included((position_2, u32::MAX))))
                .map(|(_, &node_id)| node_id)
                .collect();
            let node_ids_set = self.node_ids_by_variant_id.entry(variant_id).or_insert_with(HashSet::new);
            for &ref_node_id in &ref_node_ids {
                node_ids_set.insert(ref_node_id);
            }
            let op_node_ids_set = self.node_ids_by_operation.entry(gov.clone()).or_insert_with(HashSet::new);
            for &ref_node_id in &ref_node_ids {
                op_node_ids_set.insert(ref_node_id);
            }
        }

        //  Do not include other operation types

        None
    }

    pub fn annotate_reference_nodes(
        &mut self,
        chromosome: &str,
        start: u32,
        end: u32,
        key: &str,
        value: Arc<dyn Any + Send + Sync + 'static>
    ) {
        // Step 1. Split reference nodes at start and end positions
        let mut split_reference_at_position = |chromosome: &str, position: u32| -> usize {
            let ref_node_id: usize = self.query_reference_node(chromosome, position)
                .unwrap_or_else(|| {
                    panic!(
                        "Expected to have a reference node containing the position {}:{}",
                        chromosome, position
                    )
                });
            self.split_reference_node(ref_node_id, position)
        };
        let start_ref_node_id: usize = split_reference_at_position(chromosome, start);
        let end_ref_node_id: usize = split_reference_at_position(chromosome, end);

        // Step 2. Add annotations to reference nodes
        let ref_node_ids: Vec<usize> = self.reference_nodes_index
            .get(chromosome)
            .unwrap()
            .range((Included((start, 0u32)), Included((end, u32::MAX))))
            .map(|(_, &node_id)| node_id)
            .collect();
        for ref_node_id in ref_node_ids {
            self.get_node_data_mut(ref_node_id).as_reference_mut().add_annotation(key, value.clone());
        }
    }

    pub fn find_genome_paths(&self, target_node_ids: &HashSet<usize>, excluded_node_ids: &HashSet<usize>) -> HashSet<VarGraphPath> {
        log_info!("Started finding paths for {} node IDs.", target_node_ids.len());

        // Step 1. Collapse the variation graph
        log_info!("\tCollapsing the variation graph.");
        let (vargraph, cyclic_paths) = self.collapse(target_node_ids);

        // Step 2. Get all node IDs in the collapsed graph
        let collapsed_node_ids: &HashSet<usize> = vargraph.get_node_ids();

        // Step 3. Identify intersecting node IDs
        let mut isec_node_ids: HashSet<usize> = target_node_ids.intersection(collapsed_node_ids).cloned().collect();

        // Step 4. Extend paths
        log_info!("\tExtending paths");
        let mut terminal_states: Vec<VarGraphTraversalState> = Vec::new();
        while isec_node_ids.len() > 0 {
            let node_id: usize = isec_node_ids.iter().next().unwrap().clone();
            log_info!("\t\tExtending paths for node ID: {}. Remaining number of nodes: {}.", node_id, isec_node_ids.len());
            let mut queue: Vec<VarGraphTraversalState> = Vec::new();
            queue.push(VarGraphTraversalState {
                visited: HashSet::from([node_id]),
                path: VecDeque::from([node_id]),
                prev_edge: None
            });

            // Extend in the forward (downstream) direction
            log_info!("\t\tExtending in the downstream direction.");
            let mut states_forward: Vec<VarGraphTraversalState> = vargraph.extend_paths(
                queue,
                target_node_ids,
                excluded_node_ids,
                VarGraphOrientations::Downstream,
                None
            );

            // Extend in the backward (upstream) direction
            log_info!("\t\tExtending in the upstream direction.");
            for state in states_forward.iter_mut() {
                state.prev_edge = None;
            }
            let states_backward: Vec<VarGraphTraversalState> = vargraph.extend_paths(
                states_forward,
                target_node_ids,
                excluded_node_ids,
                VarGraphOrientations::Upstream,
                None
            );

            for state in states_backward {
                for node_id in state.path.iter() {
                    if cyclic_paths.contains_key(node_id) {
                        for nested_node_id in cyclic_paths.get(node_id).unwrap() {
                            isec_node_ids.remove(nested_node_id);
                        }
                    }
                    isec_node_ids.remove(node_id);
                }
                terminal_states.push(state);
            }
        }

        // Step 6. Replace cyclic variant node IDs with the node IDs pertaining to the path
        log_info!("\tReplacing cyclic node IDs with the nested node IDs.");
        for state in terminal_states.iter_mut() {
            // Replace cyclic variant node IDs with the node IDs pertaining to the path
            'replace_variant_node_loop: loop {
                let mut prev_node_id: Option<usize> = None;
                for i in 0..state.path.len() {
                    if i > 0 {
                        prev_node_id = Some(*state.path.get(i - 1).unwrap());
                    }
                    let curr_node_id: usize = *state.path.get(i).unwrap();
                    if cyclic_paths.contains_key(&curr_node_id) {
                        state.path.remove(i);
                        let tail: VecDeque<usize> = state.path.split_off(i);
                        let mut cyclic_path: Vec<usize> = cyclic_paths.get(&curr_node_id).unwrap().clone();

                        // Depending on the path direction, reverse the cyclic path
                        if let Some(prev_node_id) = prev_node_id {
                            let edge: &VarGraphEdge = vargraph.get_edge_data(prev_node_id, curr_node_id);
                            if edge.orientation_from == VarGraphOrientations::Upstream {
                                cyclic_path.reverse();
                            }
                        }

                        state.path.extend(cyclic_path);
                        state.path.extend(tail);
                        continue 'replace_variant_node_loop;
                    }
                }
                break 'replace_variant_node_loop;
            }

            // Insert cyclic variant node IDs
            'insert_variant_node_loop: loop {
                for i in 0..state.path.len() - 1 {
                    let curr_node_id: usize = state.path.get(i).unwrap().clone();
                    let next_node_id: usize = state.path.get(i + 1).unwrap().clone();
                    let curr_variant_node_ids_set: HashSet<usize> = self.get_outgoing_variant_node_ids(curr_node_id);
                    let next_variant_node_ids_set: HashSet<usize> = self.get_outgoing_variant_node_ids(next_node_id);
                    for &variant_node_id in curr_variant_node_ids_set.intersection(&next_variant_node_ids_set) {
                        if self.is_cyclic_variant_node(variant_node_id) {
                            state.path.insert(i + 1, variant_node_id);
                            continue 'insert_variant_node_loop;
                        }
                    }
                }
                break 'insert_variant_node_loop;
            }
        }

        // Step 6. Convert VarGraphTraversalState to VarGraphPath
        log_info!("\tConverting traversal states to paths.");
        let mut vargraph_paths: HashSet<VarGraphPath> = HashSet::new();
        let terminal_states_set: HashSet<VarGraphTraversalState> = terminal_states.into_iter().collect();
        for state in terminal_states_set.iter() {
            let mut path: VarGraphPath = VarGraphPath::new();
            match state.path.len() {
                0 => {
                    panic!("Expected to have at least 1 node in the path. VarGraphTraversalState: {:#?}", state);
                },
                1 => {
                    let node_id: usize = *state.path.get(0).unwrap();
                    let node: &VarGraphNode = self.get_node_data(node_id);
                    let end: u32 = node.get_sequence_length().saturating_sub(1);
                    let segment: VarGraphSegment = VarGraphSegment::new(
                        node_id,
                        0,
                        end,
                        &*node.get_sequence(),
                        None,
                        None,
                        None,
                        None
                    );
                    path.push_back(segment);
                },
                2 => {
                    let first_node_id: usize = *state.path.get(0).unwrap();
                    let second_node_id: usize = *state.path.get(1).unwrap();
                    let first_node: &VarGraphNode = self.get_node_data(first_node_id);
                    let second_node: &VarGraphNode = self.get_node_data(second_node_id);
                    let edge: &VarGraphEdge = self.get_edge_data(first_node_id, second_node_id);
                    let segment_1: VarGraphSegment = VarGraphSegment::new(
                        first_node_id,
                        0,
                        first_node.get_sequence_length().saturating_sub(1),
                        &*first_node.get_sequence(),
                        None,
                        Some(second_node_id),
                        None,
                        Some(edge.orientation_from.clone())
                    );
                    let segment_2: VarGraphSegment = VarGraphSegment::new(
                        second_node_id,
                        0,
                        second_node.get_sequence_length().saturating_sub(1),
                        &*second_node.get_sequence(),
                        Some(first_node_id),
                        None,
                        Some(edge.orientation_to.clone()),
                        None
                    );
                    path.push_back(segment_1);
                    path.push_back(segment_2);
                },
                _ => {
                    // First node
                    let first_node_id: usize = *state.path.get(0).unwrap();
                    let second_node_id: usize = *state.path.get(1).unwrap();
                    let first_node: &VarGraphNode = self.get_node_data(first_node_id);
                    let edge: &VarGraphEdge = self.get_edge_data(first_node_id, second_node_id);
                    let end: u32 = first_node.get_sequence_length().saturating_sub(1);
                    let segment: VarGraphSegment = VarGraphSegment::new(
                        first_node_id,
                        0,
                        end,
                        &*first_node.get_sequence(),
                        None,
                        Some(second_node_id),
                        None,
                        Some(edge.orientation_from.clone())
                    );
                    path.push_back(segment);

                    // Second node to the second last node
                    for i in 1..=state.path.len() - 2 {
                        let prev_node_id: usize = *state.path.get(i - 1).unwrap();
                        let curr_node_id: usize = *state.path.get(i).unwrap();
                        let next_node_id: usize = *state.path.get(i + 1).unwrap();
                        let curr_node: &VarGraphNode = self.get_node_data(curr_node_id);
                        let prev_edge: &VarGraphEdge = self.get_edge_data(prev_node_id, curr_node_id);
                        let next_edge: &VarGraphEdge = self.get_edge_data(curr_node_id, next_node_id);
                        let end: u32 = curr_node.get_sequence_length().saturating_sub(1);
                        let segment: VarGraphSegment = VarGraphSegment::new(
                            curr_node_id,
                            0,
                            end,
                            &*curr_node.get_sequence(),
                            Some(prev_node_id),
                            Some(next_node_id),
                            Some(prev_edge.orientation_to.clone()),
                            Some(next_edge.orientation_from.clone())
                        );

                        path.push_back(segment);
                    }

                    // Last node
                    let last_node_id: usize = *state.path.get(state.path.len() - 1).unwrap();
                    let second_last_node_id: usize = *state.path.get(state.path.len() - 2).unwrap();
                    let last_node: &VarGraphNode = self.get_node_data(last_node_id);
                    let edge: &VarGraphEdge = self.get_edge_data(second_last_node_id, last_node_id);
                    let end: u32 = last_node.get_sequence_length().saturating_sub(1);
                    let segment: VarGraphSegment = VarGraphSegment::new(
                        last_node_id,
                        0,
                        end,
                        &*last_node.get_sequence(),
                        Some(second_last_node_id),
                        None,
                        Some(edge.orientation_to.clone()),
                        None
                    );
                    path.push_back(segment);
                }
            }

            vargraph_paths.insert(path);
        }

        log_info!("Finished finding {} unique paths for {} node IDs.", vargraph_paths.len(), target_node_ids.len());

        vargraph_paths
    }

    pub fn find_transcript_path(
        &self,
        transcript_model_id: usize,
        reference_transcript_ids: &str,
        node_ids: &HashSet<usize>
    ) -> VarGraphPath {
        log_info!("Started finding transcript path for {} node IDs.", node_ids.len());

        // Get all exonic reference node IDs
        let mut exonic_reference_node_ids: HashSet<usize> = node_ids.clone();
        let mut non_exonic_reference_node_ids: HashSet<usize> = HashSet::new();
        for &reference_node_id in self.get_reference_node_ids().iter() {
            if let Some(reference_node_type) = self.get_reference_node_type(reference_node_id) {
                if reference_node_type == VarGraphReferenceNodeTypes::Exonic {
                    exonic_reference_node_ids.insert(reference_node_id);
                } else {
                    non_exonic_reference_node_ids.insert(reference_node_id);
                }
            } else {
                non_exonic_reference_node_ids.insert(reference_node_id);
            }
        }

        // Target node IDs are a union of the exonic reference node IDs and the node IDs provided
        let target_node_ids: HashSet<usize> = exonic_reference_node_ids.union(node_ids).cloned().collect();

        // Remove "bridge" reference node IDs from the non_exonic_reference_node_ids
        // A reference node ID is a bridge node if it is an intronic node that is specified
        // as part of an SNV/MNV/INS flanking a splicing.
        // For example:
        //  We have the following transcript structure:
        //      chrA:5:D:chrA:10:U (non-canonical splicing)
        //      chrA:10:D:CCC:chrA:11:U (insertion)
        //      chrA:11:I:chrA:15:I (reference exon)
        //  Then chrA:10 is considered a bridge reference node; it should be included in path traversal
        //  but excluded from the path sequence.
        let mut bridge_node_ids: HashSet<usize> = HashSet::new();
        let mut variant_ids: Vec<usize> = self.get_variant_ids().into_iter().collect();
        variant_ids.sort();
        for variant_id in variant_ids.iter() {
            let gov: &GraphOperationView = self.operation_by_variant_id.get(variant_id).unwrap();
            if *gov.get_operation_type_1() == GraphOperationType::Downstream &&
                *gov.get_operation_type_2() == GraphOperationType::Upstream {
                let reference_node_id_1: usize = self.query_reference_node(gov.get_chromosome_1(), gov.get_position_1()).unwrap();
                let reference_node_id_2: usize = self.query_reference_node(gov.get_chromosome_2(), gov.get_position_2()).unwrap();
                if non_exonic_reference_node_ids.contains(&reference_node_id_1) {
                    bridge_node_ids.insert(reference_node_id_1);
                }
                if non_exonic_reference_node_ids.contains(&reference_node_id_2) {
                    bridge_node_ids.insert(reference_node_id_2);
                }
            }
        }
        non_exonic_reference_node_ids.retain(|x| !bridge_node_ids.contains(x));

        let paths: HashSet<VarGraphPath> = self.find_genome_paths(
            &target_node_ids,
            &non_exonic_reference_node_ids
        );

        assert!(
            paths.len() == 1,
            "Expected to find exactly 1 path. Found {} paths for transcript model ID {} and reference transcript IDs {}.",
            paths.len(), transcript_model_id, reference_transcript_ids
        );

        let mut path: VarGraphPath = paths.into_iter().next().unwrap();

        // Disable bridge node segments in the path
        for segment in path.segments.iter_mut() {
            if bridge_node_ids.contains(&segment.node_id) {
                segment.disable();
            }
        }

        // If the first valid variant (i.e. Downstream / Upstream / Include operation)
        // does not match the first segment of the path, reverse the path.
        let mut first_enabled_segment_index: usize = 0;
        for (i, segment) in path.segments.iter().enumerate() {
            if segment.is_enabled() {
                first_enabled_segment_index = i;
                break;
            }
        }
        let mut path_was_reversed = false;
        for variant_id in variant_ids.iter() {
            let gov: &GraphOperationView = self.operation_by_variant_id.get(variant_id).unwrap();
            if *gov.get_operation_type_1() == GraphOperationType::Downstream ||
                *gov.get_operation_type_1() == GraphOperationType::Upstream ||
                *gov.get_operation_type_1() == GraphOperationType::Include {
                let node_ids: &HashSet<usize> = self.node_ids_by_variant_id.get(variant_id).unwrap();
                if node_ids.contains(&path.segments[first_enabled_segment_index].node_id) == false {
                    path.segments.make_contiguous().reverse();
                    for segment in path.segments.iter_mut() {
                        std::mem::swap(&mut segment.entry_orientation, &mut segment.exit_orientation);
                        std::mem::swap(&mut segment.prev_node_id, &mut segment.next_node_id);
                    }
                    path_was_reversed = true;
                }
                break;
            }
        }

        // Assign entry_orientation based on the transcript strand (from the 0th
        // index operation, which is the starting point of the transcript).
        // For reverse-strand transcripts, segments must be in descending genomic
        // order (high to low) with per-segment RC via entry_orientation = Downstream.
        // If the direction correction above didn't reverse the path, do it now.
        // For forward-strand transcripts, the existing orientations are already
        // correct — no override needed.
        let transcript_strand: Strand = self.operation_by_variant_id
            .get(&variant_ids[0])
            .map(|gov| gov.get_strand_1().clone())
            .unwrap_or(Strand::Forward);

        if transcript_strand == Strand::Reverse {
            if !path_was_reversed {
                path.segments.make_contiguous().reverse();
                for segment in path.segments.iter_mut() {
                    std::mem::swap(&mut segment.entry_orientation, &mut segment.exit_orientation);
                    std::mem::swap(&mut segment.prev_node_id, &mut segment.next_node_id);
                }
            }
            for segment in path.segments.iter_mut() {
                segment.entry_orientation = Some(VarGraphOrientations::Downstream);
            }
        }

        path
    }

    pub fn get_incoming_node_ids(&self, id: usize) -> &HashSet<usize> {
        self.multidigraph.get_incoming_node_ids(id)
    }

    pub fn get_outgoing_node_ids(&self, id: usize) -> &HashSet<usize> {
        self.multidigraph.get_outgoing_node_ids(id)
    }

    pub fn get_outgoing_reference_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
        for &outgoing_node_id in self.get_outgoing_node_ids(id) {
            if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Reference {
                outgoing_node_ids.insert(outgoing_node_id);
            }
        }
        outgoing_node_ids
    }

    pub fn get_outgoing_variant_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
        for &outgoing_node_id in self.get_outgoing_node_ids(id) {
            if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Variant {
                outgoing_node_ids.insert(outgoing_node_id);
            }
        }
        outgoing_node_ids
    }

    /// Displays the MultiDiGraph.
    pub fn display(&self) {
        println!("-------------------------------------------");
        println!("VarGraph ({} nodes, {} edges)", self.multidigraph.node_ids.len(), self.multidigraph.outgoing_edges.keys().len());
        println!("-------------------------------------------");
        println!("Node IDs");
        println!("-------------------------------------------");
        for &node_id in self.multidigraph.node_ids.iter() {
            let node: &VarGraphNode = self.get_node_data(node_id);
            println!("{:8}|{}|{}:{}-{}:{}",
                     node_id,
                     self.get_node_type(node_id).as_str(),
                     node.get_chromosome_1(),
                     node.get_position_1(),
                     node.get_chromosome_2(),
                     node.get_position_2()
            );
            // println!("{:8}|{}|{}:{}-{}:{}|{}",
            //          node_id,
            //          self.get_node_type(*node_id).as_str(),
            //          vargraph_node.get_chromosome_1(),
            //          vargraph_node.get_position_1(),
            //          vargraph_node.get_chromosome_2(),
            //          vargraph_node.get_position_2(),
            //          vargraph_node.get_sequence()
            // );
        }
        println!("-----------------------------------------------------------");
        println!("From                 | To                  | Enabled       |");
        println!("-----------------------------------------------------------");
        for (&from_node_id,destinations) in self.multidigraph.outgoing_edges.iter() {
            for &to_node_id in destinations.iter() {
                let edge: &VarGraphEdge = self.get_edge_data(from_node_id, to_node_id);
                println!("{:20} | {:20} | {:10}", from_node_id, to_node_id, edge.enabled);
            }
        }
        println!("-----------------------------------------------------------");
    }

    /// Get the total number of nodes.
    ///
    /// Returns:
    /// - Total number of nodes.
    pub fn get_nodes_count(&self) -> usize {
        self.multidigraph.node_ids.len()
    }

    pub fn get_node_ids(&self) -> &HashSet<usize> {
        &self.multidigraph.node_ids
    }

    /// Get reference node Ids.
    pub fn get_reference_node_ids(&self) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        for (chromosome, tree) in self.reference_nodes_index.iter() {
            for (key, value) in tree.iter() {
                node_ids.push(*value);
            }
        }
        node_ids
    }

    pub fn get_source_reference_node_ids(&self) -> Vec<usize> {
        let mut source_node_ids: Vec<usize> = Vec::new();
        for &node_id in self.get_node_ids() {
            if self.get_node_type(node_id) == VarGraphNodeTypes::Reference {
                let outgoing_node_ids: &HashSet<usize> = self.get_outgoing_node_ids(node_id);
                let mut has_outgoing_upstream: bool = false;
                let mut has_outgoing_downstream: bool = false;
                for &outgoing_node_id in outgoing_node_ids {
                    let edge: &VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id);
                    if edge.orientation_from == VarGraphOrientations::Upstream {
                        has_outgoing_upstream = true;
                    }
                    if edge.orientation_from == VarGraphOrientations::Downstream {
                        has_outgoing_downstream = true;
                    }
                }
                if has_outgoing_downstream == true && has_outgoing_upstream == false {
                    source_node_ids.push(node_id);
                }
            }
        }
        source_node_ids
    }

    pub fn get_sink_node_ids(&self) -> Vec<usize> {
        let mut sink_node_ids: Vec<usize> = Vec::new();
        for &node_id in self.get_node_ids().iter() {
            let outgoing_node_ids = self.get_outgoing_node_ids(node_id);
            let mut has_outgoing_upstream: bool = false;
            let mut has_outgoing_downstream: bool = false;
            for &outgoing_node_id in outgoing_node_ids {
                let edge: &VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id);
                if edge.orientation_from == VarGraphOrientations::Upstream {
                    has_outgoing_upstream = true;
                }
                if edge.orientation_from == VarGraphOrientations::Downstream {
                    has_outgoing_downstream = true;
                }
            }
            if has_outgoing_downstream == false && has_outgoing_upstream == true {
                sink_node_ids.push(node_id);
            }
        }
        sink_node_ids
    }

    pub fn get_variant_ids(&self) -> HashSet<usize> {
        self.operation_by_variant_id.keys().cloned().collect()
    }

    pub fn get_variant_node_ids(&self) -> Vec<usize> {
        let mut variant_node_ids: Vec<usize> = Vec::new();
        for (&variant_id, node_ids) in self.node_ids_by_variant_id.iter() {
            for &node_id in node_ids.iter() {
                if self.get_node_type(node_id) == VarGraphNodeTypes::Variant {
                    variant_node_ids.push(node_id);
                }
            }
        }
        variant_node_ids.sort();
        variant_node_ids
    }

    pub fn merge(vargraphs: Vec<VarGraph>) -> VarGraph {
        log_info!("Started merging {} VarGraph objects into one.", vargraphs.len());

        let mut merged_vargraph: VarGraph = VarGraph::new();

        // Step 1. Assign a new node ID for each node
        let mut new_node_ids_map: HashMap<(usize, usize), usize> = HashMap::new();
        for (i, vargraph) in vargraphs.iter().enumerate() {
            let mut sorted: Vec<usize> = vargraph.multidigraph.node_ids.iter().cloned().collect();
            sorted.sort();
            for node_id in sorted {
                new_node_ids_map.insert((i, node_id), merged_vargraph.multidigraph.next_node_id);
                merged_vargraph.multidigraph.next_node_id += 1;
            }
        }

        // Step 2. Merge the MultiDiGraph objects
        let mut vargraph_id: usize = 0;
        for vargraph in vargraphs {
            // multigraph.node_ids
            for &node_id in vargraph.multidigraph.node_ids.iter() {
                let new_node_id: usize = *new_node_ids_map.get(&(vargraph_id, node_id)).unwrap();
                merged_vargraph.multidigraph.node_ids.insert(new_node_id);
            }

            // multigraph.outgoing_edges
            for (&source_node_id, outgoing_node_ids) in vargraph.multidigraph.outgoing_edges.iter() {
                let new_source_node_id: usize = *new_node_ids_map.get(&(vargraph_id, source_node_id)).unwrap();
                let entry = merged_vargraph
                    .multidigraph
                    .outgoing_edges
                    .entry(new_source_node_id)
                    .or_default();
                for &outgoing_node_id in outgoing_node_ids.iter() {
                    let new_outgoing_node_id: usize = *new_node_ids_map.get(&(vargraph_id, outgoing_node_id)).unwrap();
                    entry.insert(new_outgoing_node_id);
                }
            }

            // multigraph.incoming_edges
            for (&target_node_id, incoming_node_ids) in vargraph.multidigraph.incoming_edges.iter() {
                let new_target_node_id: usize = *new_node_ids_map.get(&(vargraph_id, target_node_id)).unwrap();
                let entry = merged_vargraph
                    .multidigraph
                    .incoming_edges
                    .entry(new_target_node_id)
                    .or_default();
                for &incoming_node_id in incoming_node_ids.iter() {
                    let new_incoming_node_id: usize = *new_node_ids_map.get(&(vargraph_id, incoming_node_id)).unwrap();
                    entry.insert(new_incoming_node_id);
                }
            }

            // multigraph.node_data
            for (node_id, data) in vargraph.multidigraph.node_data {
                let new_node_id: usize = *new_node_ids_map.get(&(vargraph_id, node_id)).unwrap();
                merged_vargraph.multidigraph.node_data.insert(new_node_id, data);
            }

            // multigraph.edge_data
            for ((source_node_id, destination_node_id), data) in vargraph.multidigraph.edge_data {
                let new_source_node_id: usize = *new_node_ids_map.get(&(vargraph_id, source_node_id)).unwrap();
                let new_destination_node_id: usize = *new_node_ids_map.get(&(vargraph_id, destination_node_id)).unwrap();
                merged_vargraph.multidigraph.edge_data.insert((new_source_node_id, new_destination_node_id), data);
            }

            // vargraph.reference_nodes_index
            for (chromosome, tree) in vargraph.reference_nodes_index {
                for ((start, end), node_id) in tree {
                    let new_node_id: usize = *new_node_ids_map.get(&(vargraph_id, node_id)).unwrap();
                    merged_vargraph.reference_nodes_index
                        .entry(chromosome.clone())
                        .or_insert(BTreeMap::new())
                        .insert((start, end), new_node_id);
                }
            }

            // vargraph.variants
            for (variant_id, gov) in vargraph.operation_by_variant_id {
                merged_vargraph.operation_by_variant_id.insert(variant_id, gov);
            }

            // vargraph.variant_ids_index
            for (variant_id, node_ids) in vargraph.node_ids_by_variant_id {
                for node_id in node_ids {
                    let new_node_id: usize = *new_node_ids_map.get(&(vargraph_id, node_id)).unwrap();
                    merged_vargraph.node_ids_by_variant_id.entry(variant_id).or_insert_with(HashSet::new).insert(new_node_id);
                }
            }

            // vargraph.variant_operations_index
            for (gov, node_ids) in vargraph.node_ids_by_operation {
                for node_id in node_ids {
                    let new_node_id: usize = *new_node_ids_map.get(&(vargraph_id, node_id)).unwrap();
                    merged_vargraph.node_ids_by_operation.entry(gov.clone()).or_insert_with(HashSet::new).insert(new_node_id);
                }
            }

            vargraph_id += 1;
        }

        log_info!("Finished merging VarGraph objects.");

        merged_vargraph
    }

    pub fn stitch_adjacent_variants(&mut self) {
        // Step 1. Identify all chromosomes harboring variant IDs
        let mut variant_chromosomes_map: HashMap<Box<str>, HashSet<usize>> = HashMap::new();
        for (&variant_id, gov) in self.operation_by_variant_id.iter() {
            variant_chromosomes_map
                .entry(gov.get_chromosome_1().into())
                .or_insert_with(HashSet::new)
                .insert(variant_id);
            variant_chromosomes_map
                .entry(gov.get_chromosome_2().into())
                .or_insert_with(HashSet::new)
                .insert(variant_id);
        }

        // Step 2. Identify overlapping variant IDs using sweep line algorithm
        // Order of variant ID in pairs tuples indicates the edge direction (5' --> 3')
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        let mut variants_map: HashMap<usize, &GraphOperationView> = HashMap::new();
        for (chromosome, variant_ids) in variant_chromosomes_map.iter() {
            // Prepare segments
            let mut segments: Vec<(u32,u32,usize)> = Vec::new(); // (start, end, variant_id)
            for &variant_id in variant_ids.iter() {
                // Only variants with downstream and upstream operations are subject to variant stitching
                let gov: &GraphOperationView = self.get_graph_operation_view(variant_id);
                if (*gov.get_operation_type_1() == GraphOperationType::Downstream || *gov.get_operation_type_1() == GraphOperationType::Upstream) &&
                    (*gov.get_operation_type_2() == GraphOperationType::Downstream || *gov.get_operation_type_2() == GraphOperationType::Upstream) {
                    if gov.get_chromosome_1() == gov.get_chromosome_2() {
                        segments.push((gov.get_position_1(), gov.get_position_2(), variant_id));
                    } else {
                        if gov.get_chromosome_1() == &**chromosome {
                            segments.push((gov.get_position_1(), gov.get_position_1(), variant_id));
                        } else if gov.get_chromosome_2() == &**chromosome {
                            segments.push((gov.get_position_2(), gov.get_position_2(), variant_id));
                        } else {
                            panic!("Invalid variant indexed for chromosome: {}.", chromosome);
                        }
                    }
                    variants_map.insert(variant_id, gov);
                }
            }

            // Sort by start ascending, tie by end ascending
            segments.sort_by_key(|&(start, end, _)| (start, end));

            // Active set keyed by end -> indices of segments with that end
            let mut active: BTreeMap<u32, Vec<usize>> = BTreeMap::new();

            // Get pairs of overlapping segments - sweep line algorithm
            for (s, e, variant_id_1) in segments {
                // Remove segments with end < current start
                let expired_ends: Vec<u32> = active
                    .range(..s) // keys < s
                    .map(|(end, _)| *end)
                    .collect();

                for end in expired_ends {
                    active.remove(&end);
                }

                let gov_1: &GraphOperationView = self.get_graph_operation_view(variant_id_1);
                for (_end, variant_ids) in active.iter() {
                    for &variant_id_2 in variant_ids.iter() {
                        let gov_2: &GraphOperationView = self.get_graph_operation_view(variant_id_2);

                        // Determine how the adjacency should be directed (5' A --> B 3')
                        if gov_2.get_position_2() == gov_1.get_position_1() + 1 {
                            pairs.push((variant_id_2, variant_id_1));
                            continue;
                        }
                        if gov_1.get_position_2() == gov_2.get_position_1() + 1 {
                            pairs.push((variant_id_1, variant_id_2));
                            continue;
                        }
                    }
                }

                // Add current segment to active
                active.entry(e).or_default().push(variant_id_1);
            }
        }

        // Step 2. Add edges between variant nodes and disable edges to reference nodes
        for (variant_id_a, variant_id_b) in pairs.iter() {
            let gov_a: &GraphOperationView = self.get_graph_operation_view(*variant_id_a);
            let gov_b: &GraphOperationView = self.get_graph_operation_view(*variant_id_b);
            let variant_node_a: usize = self.get_variant_node_id(*variant_id_a);
            let variant_node_b: usize = self.get_variant_node_id(*variant_id_b);
            let outgoing_node_id_a: usize = self.query_reference_node(gov_a.get_chromosome_2(), gov_a.get_position_2()).unwrap();
            let outgoing_node_id_b: usize = self.query_reference_node(gov_b.get_chromosome_1(), gov_b.get_position_1()).unwrap();

            if !(gov_a.get_variant_type() == VariantType::Breakpoint &&
                gov_b.get_variant_type() == VariantType::Breakpoint) {
                // Add edges between variant nodes
                let edge_1: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Downstream,
                    VarGraphOrientations::Upstream
                );
                let edge_2: VarGraphEdge = VarGraphEdge::new(
                    VarGraphOrientations::Upstream,
                    VarGraphOrientations::Downstream
                );
                self.multidigraph.add_edge(variant_node_a, variant_node_b, edge_1);
                self.multidigraph.add_edge(variant_node_b, variant_node_a, edge_2);

                // Disable edges to reference nodes
                // Disable edge from Variant A to Variant A position_2 reference node
                self.disable_edge(variant_node_a, outgoing_node_id_a);

                // Disable edge from Variant B position_1 reference node to Variant B
                self.disable_edge(outgoing_node_id_b, variant_node_b);
            }
        }
    }

    // pub fn get_all_linearized_contigs(&self) -> Vec<VarGraphPath> {
    //     // Step 1. Get all linearized contigs involving all variant nodes
    //     let mut vargraph_paths: Vec<VarGraphPath> = self.get_paths(self.get_variant_node_ids());
    // 
    //     // Step 2. Identify any reference  that does not have a variant node
    //     let mut all_contig_names: HashSet<Box<str>> = HashSet::new();
    //     for reference_node_id in self.get_reference_node_ids() {
    //         let node: &dyn VarGraphNode = self.get_node(reference_node_id);
    //         all_contig_names.insert(node.get_chromosome_1().clone());
    //     }
    //     for variant_node_id in self.get_variant_node_ids().iter() {
    //         let node: &dyn VarGraphNode = self.get_node(*variant_node_id);
    //         all_contig_names.remove(&node.get_chromosome_1());
    //         all_contig_names.remove(&node.get_chromosome_2());
    //     }
    //     for contig_name in all_contig_names.iter() {
    //         let reference_node_ids: Vec<usize> = self.reference_nodes_index
    //             .iter()
    //             .filter(|((key_str, _, _), _)| &**key_str == &**contig_name)
    //             .map(|(_, value)| *value)
    //             .collect();
    //         for reference_node_id in reference_node_ids.iter() {
    //             let node: &dyn VarGraphNode = self.get_node(*reference_node_id);
    //             let mut vargraph_path: VarGraphPath = VarGraphPath::new();
    //             let vargraph_segment: VarGraphSegment = VarGraphSegment::new(
    //                 *reference_node_id,
    //                 0,
    //                 (node.get_sequence_length() as usize).saturating_sub(1),
    //                 &*node.get_sequence(),
    //                 VarGraphStrands::Forward
    //             );
    //             vargraph_path.push_back(vargraph_segment);
    //             vargraph_paths.push(vargraph_path);
    //         }
    //     }
    // 
    //     vargraph_paths
    // }
}

/// Helper methods
impl VarGraph {
    fn add_edge(&mut self, source_node_id: usize, destination_node_id: usize, data: VarGraphEdge) {
        self.multidigraph.add_edge(source_node_id, destination_node_id, data);
    }

    fn build_cyclic_paths(&self, node_ids: &HashSet<usize>) -> HashMap<usize, Vec<usize>> {
        let target_node_ids: HashSet<usize> = node_ids.clone();

        // Step 1. Build micro cyclic paths
        let mut micro_paths: HashMap<usize, Vec<usize>> = HashMap::new();
        for &node_id in node_ids.iter() {
            if self.is_cyclic_variant_node(node_id) {
                let (mut start, mut end) = (None, None);
                for &outgoing_node_id in self.get_outgoing_node_ids(node_id).iter() {
                    let edge: &VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id);
                    match edge.orientation_from {
                        VarGraphOrientations::Upstream => start = Some(outgoing_node_id),
                        VarGraphOrientations::Downstream => end = Some(outgoing_node_id),
                        _ => {}
                    }
                }
                let (start, end) = (
                    start.expect("Missing upstream (start) node"),
                    end.expect("Missing downstream (end) node")
                );
                let node: &VarGraphNode = self.get_node_data(node_id);
                let variant_node: &VarGraphVariantNode = node.as_variant();
                let num_cycles: u16 = variant_node.graph_operation_view.get_num_cycles().unwrap();
                let local_path: Vec<usize> = self.traverse_cycle_path(start, end, &target_node_ids).expect("Missing path");
                let mut local_path_cycled: Vec<usize> = Vec::new();
                for _ in 0..num_cycles {
                    local_path_cycled.extend(local_path.iter());
                }
                micro_paths.insert(node_id, local_path_cycled);
            }
        }

        // Step 2. Expand cyclic paths
        fn expand_path(
            definitions: &HashMap<usize, Vec<usize>>,
            path: Vec<usize>
        ) -> Vec<usize> {
            let mut path = path;
            for (i, node_id) in path.iter().enumerate() {
                if let Some(segment_2) = definitions.get(node_id) {
                    let mut new_path = Vec::with_capacity(path.len() - 1 + segment_2.len());
                    new_path.extend_from_slice(&path[..i]);
                    new_path.extend_from_slice(segment_2);
                    new_path.extend_from_slice(&path[i + 1..]);
                    return expand_path(definitions, new_path);
                }
            }
            path
        }

        let mut expanded_paths: HashMap<usize, Vec<usize>> = HashMap::new();
        for (&variant_id, path) in micro_paths.iter() {
            let expanded_path: Vec<usize> = expand_path(&micro_paths, path.clone());
            expanded_paths.insert(variant_id, expanded_path);
        }

        // Step 3. Consolidate paths
        fn is_path_subset(path_1: &Vec<usize>, path_2: &Vec<usize>) -> bool {
            fn is_subpath(short: &[usize], long: &[usize]) -> bool {
                if short.len() > long.len() {
                    return false;
                }
                if short.is_empty() {
                    return true;
                }
                long.windows(short.len()).any(|window| window == short)
            }

            is_subpath(path_1, path_2) || is_subpath(path_2, path_1)
        }

        fn consolidate_paths(
            paths: HashMap<usize, Vec<usize>>
        ) -> HashMap<usize, Vec<usize>> {
            let mut uf: UnionFind = UnionFind::new();
            for (&variant_id, path) in paths.iter() {
                uf.union(variant_id as u32, variant_id as u32);
            }

            for (&variant_id_1, path_1) in paths.iter() {
                for (&variant_id_2, path_2) in paths.iter() {
                    if variant_id_1 == variant_id_2 {
                        continue;
                    }
                    if is_path_subset(path_1, path_2) {
                        uf.union(variant_id_1 as u32, variant_id_2 as u32);
                    }
                }
            }

            let mut consolidated_paths: HashMap<usize, Vec<usize>> = HashMap::new();
            for cluster in uf.get_clusters() {
                let variant_node_id: usize = cluster
                    .iter()
                    .max_by_key(|id| paths.get(&(**id as usize)).unwrap().len())
                    .copied()
                    .unwrap() as usize;
                let path = paths.get(&variant_node_id).unwrap();
                consolidated_paths.insert(variant_node_id, path.clone());
            }

            consolidated_paths
        }

        let consolidated_paths: HashMap<usize, Vec<usize>> = consolidate_paths(expanded_paths);

        consolidated_paths
    }

    fn collapse(&self, node_ids: &HashSet<usize>) -> (VarGraph, HashMap<usize, Vec<usize>>) {
        // Step 1. Identify cyclic nodes
        let mut cyclic_node_ids: HashSet<usize> = HashSet::new();
        for node_id in node_ids.iter() {
            if self.is_cyclic_variant_node(*node_id) {
                cyclic_node_ids.insert(*node_id);
            }
        }

        // Step 2. Build macro paths
        let cyclic_paths: HashMap<usize, Vec<usize>> = self.build_cyclic_paths(&node_ids);

        // Step 3. Clone the VarGraph
        let mut vargraph: VarGraph = self.clone();

        // Step 4. Rewrite the variation graph with cyclic paths
        for (variant_node_id, path) in cyclic_paths.iter() {
            let path_set: HashSet<usize> = path.iter().copied().collect();
            let start: usize = *path.first().unwrap();
            for &outgoing_node_id in self.get_outgoing_node_ids(start) {
                if self.is_edge_enabled(start, outgoing_node_id) == false {
                    continue;
                }
                if outgoing_node_id != *variant_node_id &&
                    path_set.contains(&outgoing_node_id) == false {
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Upstream,
                        VarGraphOrientations::Downstream
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Downstream,
                        VarGraphOrientations::Upstream
                    );
                    vargraph.add_edge(*variant_node_id, outgoing_node_id, edge_1);
                    vargraph.add_edge(outgoing_node_id, *variant_node_id, edge_2);
                }
            }

            let end: usize = *path.last().unwrap();
            for &outgoing_node_id in self.get_outgoing_node_ids(end) {
                if self.is_edge_enabled(end, outgoing_node_id) == false {
                    continue;
                }
                if outgoing_node_id != *variant_node_id &&
                    path_set.contains(&outgoing_node_id) == false {
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Downstream,
                        VarGraphOrientations::Upstream
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        VarGraphOrientations::Upstream,
                        VarGraphOrientations::Downstream
                    );
                    vargraph.add_edge(*variant_node_id, outgoing_node_id, edge_1);
                    vargraph.add_edge(outgoing_node_id, *variant_node_id, edge_2);
                }
            }
        }

        // Step 5. Remove nodes contained in each cyclic path
        for (_, path) in cyclic_paths.iter() {
            for node_id in path.iter() {
                vargraph.remove_node(*node_id);
            }
        }

        // Step 6. Remove cyclic nodes that are collapsed
        for node_id in cyclic_node_ids.iter() {
            if cyclic_paths.contains_key(node_id) == false {
                vargraph.remove_node(*node_id);
            }
        }

        (vargraph, cyclic_paths)
    }

    fn disable_edge(&mut self, from: usize, to: usize) {
        self.get_edge_data_mut(from, to).enabled = false;
        self.get_edge_data_mut(to, from).enabled = false;
    }

    fn extend_cycle_path(
        &self,
        initial_state: &VarGraphTraversalState,
        from: usize,
        to: usize,
        target_node_ids: &HashSet<usize>,
        direction: VarGraphOrientations
    ) -> Option<VarGraphTraversalState> {
        assert_eq!(target_node_ids.contains(&to), true);

        let mut state: VarGraphTraversalState = initial_state.clone();

        'extension_loop: loop {
            let curr_node_id: usize = match direction {
                VarGraphOrientations::Downstream => *state.path.back().unwrap(),
                VarGraphOrientations::Upstream => *state.path.front().unwrap(),
                _ => panic!("Invalid direction")
            };

            // If we're already at the destination, treat this as terminal and don't extend
            if curr_node_id == to {
                return Some(state);
            }

            // Check if the current node is an entry or exit to a cycle
            let mut cyclic_variant_node_id: Option<usize> = None;
            if curr_node_id != from && curr_node_id != to {
                for &variant_node_id in self.get_outgoing_variant_node_ids(curr_node_id).iter() {
                    if self.is_cyclic_variant_node(variant_node_id) {
                        cyclic_variant_node_id = Some(variant_node_id);
                        break;
                    }
                }
            }

            // Update the state
            let mut target_states: Vec<VarGraphTraversalState> = Vec::new();
            let mut all_states: Vec<VarGraphTraversalState> = Vec::new();
            if cyclic_variant_node_id.is_some() {
                // Find the cyclic path start or end (mate) node ID
                let variant_node_id: usize = cyclic_variant_node_id.unwrap();
                let mut mate_node_id: Option<usize> = None;
                for &outgoing_node_id in self.get_outgoing_node_ids(variant_node_id).iter() {
                    if outgoing_node_id != curr_node_id {
                        mate_node_id = Some(outgoing_node_id);
                        break;
                    }
                }
                assert!(mate_node_id.is_some());
                let mate_node_id: usize = mate_node_id.unwrap();

                for &outgoing_node_id in self.get_outgoing_node_ids(mate_node_id).iter() {
                    if self.is_edge_enabled(mate_node_id, outgoing_node_id) == false {
                        continue;
                    }

                    // Avoid revisiting nodes
                    if state.visited.contains(&outgoing_node_id) {
                        continue;
                    }

                    // Avoid cyclic variant nodes
                    if self.is_cyclic_variant_node(outgoing_node_id) {
                        continue;
                    }

                    // Avoid variant nodes that are not target nodes
                    if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Variant &&
                        target_node_ids.contains(&outgoing_node_id) == false {
                        continue;
                    }

                    let edge: &VarGraphEdge = self.get_edge_data(mate_node_id, outgoing_node_id);

                    // Orientation constraints
                    if state.prev_edge.is_none() {
                        if direction != edge.orientation_from {
                            continue;
                        }
                    } else {
                        if state.prev_edge.as_ref().unwrap().orientation_to == edge.orientation_from {
                            continue;
                        }
                    }

                    // Create a new traversal state
                    let mut new_state: VarGraphTraversalState = state.clone();
                    if direction == VarGraphOrientations::Downstream {
                        new_state.path.pop_back();
                        new_state.path.push_back(variant_node_id);
                        new_state.path.push_back(outgoing_node_id);
                    } else {
                        new_state.path.pop_front();
                        new_state.path.push_front(variant_node_id);
                        new_state.path.push_front(outgoing_node_id);
                    }
                    new_state.visited.insert(curr_node_id);
                    new_state.visited.insert(variant_node_id);
                    new_state.visited.insert(mate_node_id);
                    new_state.visited.insert(outgoing_node_id);
                    new_state.prev_edge = Some(edge.clone());

                    if target_node_ids.contains(&outgoing_node_id) {
                        target_states.push(new_state);
                    } else {
                        all_states.push(new_state);
                    }
                }

            } else {
                for &outgoing_node_id in self.get_outgoing_node_ids(curr_node_id) {
                    if self.is_edge_enabled(curr_node_id, outgoing_node_id) == false {
                        continue;
                    }

                    // Avoid revisiting nodes
                    if state.visited.contains(&outgoing_node_id) {
                        continue;
                    }

                    // Avoid cyclic variant nodes
                    if self.is_cyclic_variant_node(outgoing_node_id) {
                        continue;
                    }

                    // Avoid variant nodes that are not target nodes
                    if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Variant &&
                        target_node_ids.contains(&outgoing_node_id) == false {
                        continue;
                    }

                    let edge: &VarGraphEdge = self.get_edge_data(curr_node_id, outgoing_node_id);

                    // Orientation constraints
                    if state.prev_edge.is_none() {
                        if direction != edge.orientation_from {
                            continue;
                        }
                    } else {
                        if state.prev_edge.as_ref().unwrap().orientation_to == edge.orientation_from {
                            continue;
                        }
                    }

                    // Create a new traversal state
                    let mut new_state: VarGraphTraversalState = state.clone();
                    if direction == VarGraphOrientations::Downstream {
                        new_state.path.push_back(outgoing_node_id);
                    } else {
                        new_state.path.push_front(outgoing_node_id);
                    }
                    new_state.visited.insert(outgoing_node_id);
                    new_state.prev_edge = Some(edge.clone());

                    if target_node_ids.contains(&outgoing_node_id) {
                        target_states.push(new_state);
                    } else {
                        all_states.push(new_state);
                    }
                }
            }

            assert!(target_states.len() <= 1);
            assert!(all_states.len() <= 1);

            let next_states: Vec<VarGraphTraversalState> = if target_states.is_empty() {
                all_states
            } else {
                target_states
            };

            if next_states.is_empty() {
                break 'extension_loop;
            } else {
                state = next_states.get(0).unwrap().clone();
            }
        }

        None
    }

    fn extend_paths(
        &self,
        mut queue: Vec<VarGraphTraversalState>,
        target_node_ids: &HashSet<usize>,
        excluded_node_ids: &HashSet<usize>,
        direction: VarGraphOrientations,
        to: Option<usize>
    ) -> Vec<VarGraphTraversalState> {
        let mut terminal_states: Vec<VarGraphTraversalState> = Vec::new();

        if to.is_some() {
            assert_eq!(target_node_ids.contains(&to.unwrap()), true);
        }

        while let Some(state) = queue.pop() {
            let curr_node_id: usize = if direction == VarGraphOrientations::Downstream {
                *state.path.back().unwrap()
            } else {
                *state.path.front().unwrap()
            };

            // If we're already at the destination, treat this as terminal and don't extend
            if let Some(dest) = to {
                if curr_node_id == dest {
                    terminal_states.push(state);
                    continue;
                }
            }

            // Push new states to the queue
            let mut all_states: Vec<VarGraphTraversalState> = Vec::new();
            let mut traversable_states: Vec<VarGraphTraversalState> = Vec::new();
            for &outgoing_node_id in self.get_outgoing_node_ids(curr_node_id) {
                if self.is_edge_enabled(curr_node_id, outgoing_node_id) == false {
                    continue;
                }

                if excluded_node_ids.contains(&outgoing_node_id) {
                    continue;
                }

                // Avoid revisiting nodes
                if state.visited.contains(&outgoing_node_id) {
                    continue;
                }

                if self.get_node_type(outgoing_node_id) == VarGraphNodeTypes::Variant {
                    if target_node_ids.contains(&outgoing_node_id) == false {
                        continue;
                    }
                }

                let edge: &VarGraphEdge = self.get_edge_data(curr_node_id, outgoing_node_id);

                // Orientation constraints
                let prev_matches: bool = match &state.prev_edge {
                    None => direction == edge.orientation_from,
                    Some(prev_edge) => prev_edge.orientation_to != edge.orientation_from,
                };
                if !prev_matches {
                    continue;
                }

                // Create a new traversal state
                let mut new_state: VarGraphTraversalState = state.clone();
                if direction == VarGraphOrientations::Downstream {
                    new_state.path.push_back(outgoing_node_id);
                } else {
                    new_state.path.push_front(outgoing_node_id);
                }
                new_state.visited.insert(outgoing_node_id);
                new_state.prev_edge = Some(edge.clone());

                if target_node_ids.contains(&outgoing_node_id) {
                    traversable_states.push(new_state.clone());
                }
                all_states.push(new_state);
            }

            let next_states: Vec<VarGraphTraversalState> = if traversable_states.is_empty() {
                all_states
            } else {
                traversable_states
            };

            if next_states.is_empty() {
                if let Some(dest) = to {
                    if state.path.contains(&dest) {
                        terminal_states.push(state);
                    }
                } else {
                    terminal_states.push(state);
                }
            } else {
                for s in next_states {
                    queue.push(s);
                }
            }
        }

        terminal_states
    }

    fn get_edge_data(&self, from: usize, to: usize) -> &VarGraphEdge {
        self.multidigraph.get_edge_data(from, to)
    }

    fn get_edge_data_mut(&mut self, from: usize, to: usize) -> &mut VarGraphEdge {
        self.multidigraph.get_edge_data_mut(from, to)
    }

    fn get_node_data(&self, node_id: usize) -> &VarGraphNode {
        self.multidigraph.get_node_data(node_id)
    }

    fn get_node_data_mut(&mut self, node_id: usize) -> &mut VarGraphNode {
        self.multidigraph.get_node_data_mut(node_id)
    }

    fn get_node_type(&self, node_id: usize) -> VarGraphNodeTypes {
        self.get_node_data(node_id).get_type()
    }

    fn get_graph_operation_view(&self, variant_id: usize) -> &GraphOperationView {
       self.operation_by_variant_id.get(&variant_id).unwrap()
    }

    fn get_reference_node_type(&self, node_id: usize) -> Option<VarGraphReferenceNodeTypes> {
        if self.get_node_type(node_id) != VarGraphNodeTypes::Reference {
            return None;
        }
        let node: &VarGraphReferenceNode = self.get_node_data(node_id).as_reference();
        if node.annotation_key_exists(VarGraphReferenceNodeAnnotationKeys::Type.as_str()) {
            if let value = node.get_annotation(VarGraphReferenceNodeAnnotationKeys::Type.as_str()) {
                if let Some(s) = value.downcast_ref::<Box<str>>() {
                    Some(VarGraphReferenceNodeTypes::from_str(s).unwrap())
                } else {
                    panic!("Box<str> type value was expected for node ID {}", node_id);
                }
            } else {
                panic!("Annotation value 'type' was expected for node ID {}", node_id);
            }
        } else {
            Some(VarGraphReferenceNodeTypes::Intergenic)
        }
    }

    fn get_variant_node_id(&self, variant_id: usize) -> usize {
        let mut variant_node_id: HashSet<usize> = HashSet::new();
        for &node_id in self.node_ids_by_variant_id.get(&variant_id).unwrap() {
            if self.get_node_type(node_id) == VarGraphNodeTypes::Variant {
                return node_id;
            }
        }
        assert!(variant_node_id.len() == 1);
        variant_node_id.iter().copied().next().unwrap()
    }

    fn is_cyclic_variant_node(&self, node_id: usize) -> bool {
        if self.get_node_type(node_id) == VarGraphNodeTypes::Variant {
            let variant_node: &VarGraphVariantNode = self.get_node_data(node_id).as_variant();
            variant_node.graph_operation_view.get_num_cycles().is_some()
        } else {
            false
        }
    }

    fn is_edge_enabled(&self, from: usize, to: usize) -> bool {
        self.get_edge_data(from, to).enabled
    }

    fn node_exists(&self, node_id: usize) -> bool {
        self.multidigraph.node_exists(node_id)
    }

    fn query_reference_node(&self, chromosome: &str, position: u32) -> Option<usize> {
        let tree = self.reference_nodes_index.get(chromosome).unwrap();
        let lower = (0u32, 0u32);
        let upper = (position, u32::MAX);
        if let Some((&(start, end), &node_id)) = tree.range((Included(lower), Included(upper))).next_back() {
            if start <= position && position <= end {
                return Some(node_id);
            }
        }
        None
    }

    fn query_variant_node(&self, gov: &GraphOperationView) -> Option<usize> {
        if self.node_ids_by_operation.contains_key(gov) {
            let mut variant_node_ids: HashSet<usize> = HashSet::new();
            for &node_id in self.node_ids_by_operation.get(gov).unwrap() {
                if self.get_node_type(node_id) == VarGraphNodeTypes::Variant {
                    variant_node_ids.insert(node_id);
                }
            }
            assert!(variant_node_ids.len() <= 1);
            if variant_node_ids.is_empty() {
                None
            } else {
                Some(variant_node_ids.iter().copied().next().unwrap())
            }
        } else {
            None
        }
    }

    /// Removes a node.
    ///
    /// Parameters:
    /// - id    :   Node ID.
    fn remove_node(&mut self, node_id: usize) {
        if self.node_exists(node_id) == false {
            return;
        }

        if self.get_node_type(node_id) == VarGraphNodeTypes::Reference {
            let (chromosome, start, end) = {
                let reference_node: &VarGraphReferenceNode = self.get_node_data(node_id).as_reference();
                (reference_node.get_chromosome().to_string().into_boxed_str(), reference_node.get_start(), reference_node.get_end())
            };
            self.reference_nodes_index
                .get_mut(&*chromosome)
                .expect("chromosome not found in self.reference_nodes_index")
                .remove(&(start, end));
        } else {
            let variant_node: &VarGraphVariantNode = self.get_node_data(node_id).as_variant();
            self.node_ids_by_operation.remove(&variant_node.graph_operation_view.clone());
            for (_, node_ids) in self.node_ids_by_variant_id.iter_mut() {
                if node_ids.contains(&node_id) {
                    node_ids.remove(&node_id);
                }
            }
            for (_, node_ids) in self.node_ids_by_operation.iter_mut() {
                if node_ids.contains(&node_id) {
                    node_ids.remove(&node_id);
                }
            }
        }
        self.multidigraph.remove_node(node_id);
    }

    /// Splits a reference node. The result is that there will be a reference node at
    /// position 'position' with sequence length of 1.
    ///
    /// Parameters:
    /// - id        :   Reference node ID.
    /// - position  :   Position to split.
    ///
    /// Returns:
    /// Node ID at position 'position'.
    fn split_reference_node(&mut self, node_id: usize, position: u32) -> usize {
        assert!(self.get_node_type(node_id) == VarGraphNodeTypes::Reference, "Node with ID {} is not a reference node", node_id);
        assert!(self.node_exists(node_id) == true, "Node with ID {} does not exist.", node_id);

        // Step 1. Fetch the reference node object
        let reference_node: &VarGraphReferenceNode = self.get_node_data(node_id).as_reference();
        if reference_node.get_sequence().len() == 1 {
            // The node is already independent
            return node_id;
        }

        // Step 2. Get all outgoing node IDs to 5' and 3' ends from the reference node
        let outgoing_node_ids: HashSet<usize> = self.get_outgoing_node_ids(node_id).clone();

        // Step 3. Split the reference node
        let returning_node_id: usize;
        if position == reference_node.get_start() {
            // Create 2 new VarGraphReferenceNode objects
            // At position 'position'
            let mut new_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                position,
                position,
                &reference_node.get_sequence()[0..1]
            );
            new_node_1.set_annotations(reference_node.get_annotations());
            let chromosome_1: Box<str> = new_node_1.get_chromosome().into();
            let start_1: u32 = new_node_1.get_start();
            let end_1: u32 = new_node_1.get_end();

            // From (position+1) to the end of the reference node
            let mut new_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                position + 1,
                reference_node.get_end(),
                &reference_node.get_sequence()[1..]
            );
            new_node_2.set_annotations(reference_node.get_annotations());
            let chromosome_2: Box<str> = new_node_2.get_chromosome().into();
            let start_2: u32 = new_node_2.get_start();
            let end_2: u32 = new_node_2.get_end();

            // Add new nodes
            let new_node_id_1: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_1));
            let new_node_id_2: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_2));

            // Add the new nodes to the index
            let key_1: (u32, u32) = (start_1, end_1);
            self.reference_nodes_index
                .entry(chromosome_1)
                .or_insert_with(BTreeMap::new)
                .insert(key_1, new_node_id_1);
            let key_2: (u32, u32) = (start_2, end_2);
            self.reference_nodes_index
                .entry(chromosome_2)
                .or_insert_with(BTreeMap::new)
                .insert(key_2, new_node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from new_node_id_1
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Upstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_1, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_1, edge_2);
                }
            }

            // Add edges between the new nodes
            let edge_1: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Downstream,
                VarGraphOrientations::Upstream
            );
            let edge_2: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Upstream,
                VarGraphOrientations::Downstream
            );
            self.add_edge(new_node_id_1, new_node_id_2, edge_1);
            self.add_edge(new_node_id_2, new_node_id_1, edge_2);

            // Add edges from all original outgoing nodes towards the 3' end from new_node_id_2
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Downstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_2, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_2, edge_2);
                }
            }

            // Remove the old node
            self.remove_node(node_id);

            // Assign the returning node ID
            returning_node_id = new_node_id_1;
        } else if position == reference_node.get_end() {
            // Create 2 new VarGraphReferenceNode objects
            // From start of the reference node to (position-1)
            let mut new_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                reference_node.get_start(),
                position - 1,
                &reference_node.get_sequence()[0..reference_node.get_sequence().len() - 1]
            );
            new_node_1.set_annotations(reference_node.get_annotations());
            let chromosome_1: Box<str> = new_node_1.get_chromosome().into();
            let start_1: u32 = new_node_1.get_start();
            let end_1: u32 = new_node_1.get_end();

            // At position 'position'
            let mut new_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                position,
                position,
                &reference_node.get_sequence()[reference_node.get_sequence().len() - 1..]
            );
            new_node_2.set_annotations(reference_node.get_annotations());
            let chromosome_2: Box<str> = new_node_2.get_chromosome().into();
            let start_2: u32 = new_node_2.get_start();
            let end_2: u32 = new_node_2.get_end();

            // Add new nodes
            let new_node_id_1: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_1));
            let new_node_id_2: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_2));

            // Add the new nodes to the index
            let key_1: (u32, u32) = (start_1, end_1);
            self.reference_nodes_index
                .entry(chromosome_1)
                .or_insert_with(BTreeMap::new)
                .insert(key_1, new_node_id_1);
            let key_2: (u32, u32) = (start_2, end_2);
            self.reference_nodes_index
                .entry(chromosome_2)
                .or_insert_with(BTreeMap::new)
                .insert(key_2, new_node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from new_node_id_1
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Upstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_1, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_1, edge_2);
                }
            }

            // Add edges between the new nodes
            let edge_1: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Downstream,
                VarGraphOrientations::Upstream
            );
            let edge_2: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Upstream,
                VarGraphOrientations::Downstream
            );
            self.add_edge(new_node_id_1, new_node_id_2, edge_1);
            self.add_edge(new_node_id_2, new_node_id_1, edge_2);

            // Add edges from all original outgoing nodes towards the 3' end from new_node_id_2
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Downstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_2, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_2, edge_2);
                }
            }

            // Remove the old node
            self.remove_node(node_id);

            // Assign the returning node ID
            returning_node_id = new_node_id_2;
        } else {
            // Create 3 new VarGraphReferenceNode objects
            // From the start of the reference node to (position-1)
            let mut new_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                reference_node.get_start(),
                position - 1,
                &reference_node.get_sequence()[..(position - reference_node.get_start()) as usize]
            );
            new_node_1.set_annotations(reference_node.get_annotations());
            let chromosome_1: Box<str> = new_node_1.get_chromosome().into();
            let start_1: u32 = new_node_1.get_start();
            let end_1: u32 = new_node_1.get_end();

            // At position 'position'
            let mut new_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                position,
                position,
                &reference_node.get_sequence()[(position - reference_node.get_start()) as usize..(position - reference_node.get_start() + 1) as usize]
            );
            new_node_2.set_annotations(reference_node.get_annotations());
            let chromosome_2: Box<str> = new_node_2.get_chromosome().into();
            let start_2: u32 = new_node_2.get_start();
            let end_2: u32 = new_node_2.get_end();

            // From (position+1) to the end of the reference node
            let mut new_node_3: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.get_chromosome(),
                position + 1,
                reference_node.get_end(),
                &reference_node.get_sequence()[(position - reference_node.get_start() + 1) as usize..]
            );
            new_node_3.set_annotations(reference_node.get_annotations());
            let chromosome_3: Box<str> = new_node_3.get_chromosome().into();
            let start_3: u32 = new_node_3.get_start();
            let end_3: u32 = new_node_3.get_end();

            // Add new nodes
            let new_node_id_1: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_1));
            let new_node_id_2: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_2));
            let new_node_id_3: usize = self.multidigraph.add_node(VarGraphNode::Reference(new_node_3));

            // Add the new nodes to the index
            let key_1: (u32, u32) = (start_1, end_1);
            self.reference_nodes_index
                .entry(chromosome_1)
                .or_insert_with(BTreeMap::new)
                .insert(key_1, new_node_id_1);

            let key_2: (u32, u32) = (start_2, end_2);
            self.reference_nodes_index
                .entry(chromosome_2)
                .or_insert_with(BTreeMap::new)
                .insert(key_2, new_node_id_2);

            let key_3: (u32, u32) = (start_3, end_3);
            self.reference_nodes_index
                .entry(chromosome_3)
                .or_insert_with(BTreeMap::new)
                .insert(key_3, new_node_id_3);

            // Add edges from all original outgoing nodes towards the 5' end from new_node_id_1
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Upstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_1, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_1, edge_2);
                }
            }

            // Add edges between the new nodes
            let edge_1: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Downstream,
                VarGraphOrientations::Upstream
            );
            let edge_2: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Upstream,
                VarGraphOrientations::Downstream
            );
            self.add_edge(new_node_id_1, new_node_id_2, edge_1);
            self.add_edge(new_node_id_2, new_node_id_1, edge_2);

            let edge_3: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Downstream,
                VarGraphOrientations::Upstream
            );
            let edge_4: VarGraphEdge = VarGraphEdge::new(
                VarGraphOrientations::Upstream,
                VarGraphOrientations::Downstream
            );
            self.add_edge(new_node_id_2, new_node_id_3, edge_3);
            self.add_edge(new_node_id_3, new_node_id_2, edge_4);

            // Add edges from all original outgoing nodes towards the 3' end from new_node_id_3
            for &outgoing_node_id in outgoing_node_ids.iter() {
                let edge: VarGraphEdge = self.get_edge_data(node_id, outgoing_node_id).clone();
                if edge.orientation_from == VarGraphOrientations::Downstream {
                    // Add edges
                    let edge_1: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_from.clone(),
                        edge.orientation_to.clone()
                    );
                    let edge_2: VarGraphEdge = VarGraphEdge::new(
                        edge.orientation_to.clone(),
                        edge.orientation_from.clone()
                    );
                    self.add_edge(new_node_id_3, outgoing_node_id, edge_1);
                    self.add_edge(outgoing_node_id, new_node_id_3, edge_2);
                }
            }

            // Remove the old node
            self.remove_node(node_id);

            // Assign the returning node ID
            returning_node_id = new_node_id_2;
        }

        returning_node_id
    }

    fn traverse_cycle_path(
        &self,
        from: usize,
        to: usize,
        target_node_ids: &HashSet<usize>
    ) -> Option<Vec<usize>> {
        let mut target_node_ids: HashSet<usize> = target_node_ids.clone();
        target_node_ids.insert(to);

        let initial_state: VarGraphTraversalState = VarGraphTraversalState {
            visited: HashSet::from([from]),
            path: VecDeque::from([from]),
            prev_edge: None
        };

        // Extend in the forward (downstream) direction
        let state: Option<VarGraphTraversalState> = self.extend_cycle_path(
            &initial_state,
            from,
            to,
            &target_node_ids,
            VarGraphOrientations::Downstream
        );

        if state.is_some() {
            Some(state.unwrap().path.into_iter().collect())
        } else {
            let initial_state: VarGraphTraversalState = VarGraphTraversalState {
                visited: HashSet::from([from]),
                path: VecDeque::from([from]),
                prev_edge: None
            };

            // Extend in the backward (upstream) direction
            let state: Option<VarGraphTraversalState> = self.extend_cycle_path(
                &initial_state,
                from,
                to,
                &target_node_ids, VarGraphOrientations::Upstream
            );

            if state.is_some() {
                Some(state.unwrap().path.into_iter().collect())
            } else {
                None
            }
        }
    }
}

impl Clone for VarGraph {
    fn clone(&self) -> Self {
        VarGraph {
            multidigraph: self.multidigraph.clone(),
            reference_nodes_index: self.reference_nodes_index.clone(),
            node_ids_by_operation: self.node_ids_by_operation.clone(),
            node_ids_by_variant_id: self.node_ids_by_variant_id.clone(),
            operation_by_variant_id: self.operation_by_variant_id.clone()
        }
    }
}
