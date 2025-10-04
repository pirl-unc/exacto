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
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::ops::Bound::Included;

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraph {
    pub multidigraph: MultiDiGraph,

    /// Nodes in the variation graph.
    ///
    /// HashMap
    ///     Key     :   Node ID.
    ///     Value   :   SequenceGraphNode object.
    pub nodes: HashMap<usize, Box<dyn VarGraphNode>>,

    /// Edges in the variation graph.
    ///
    /// HashMap
    ///     Key     :   (node ID, node ID).
    ///     Value   :   VarGraphEdge object.
    pub edges: HashMap<(usize, usize), VarGraphEdge>,

    /// Reference nodes index.
    ///
    /// BTreeMap:
    ///     Key     :   (chromosome, start, end).
    ///     Value   :   Node ID.
    pub reference_nodes_index: BTreeMap<(Box<str>, usize, usize), usize>,

    /// Variant nodes index.
    ///
    /// BTreeMap:
    ///     Key     :   (chromosome_1, position_1, orientation_1, chromosome_2, position_2, orientation_2).
    ///     Value   :   Node ID.
    pub variant_nodes_index: BTreeMap<(Box<str>, usize, VarGraphOrientations, Box<str>, usize, VarGraphOrientations), usize>,

    /// Region annotations.
    ///
    /// Outer HashMap
    ///     Key     : (chromosome, start, end).
    ///     Value   : An inner HashMap representing the region gene_annotation.
    ///
    /// Inner HashMap
    ///     Key     : Attribute name.
    ///     Value   : Attribute value.
    pub annotations: HashMap<(Box<str>, usize, usize), HashMap<Box<str>, Box<dyn Any>>>
}

// Public methods
impl VarGraph {
    /// Constructors
    pub fn from_fasta_file(
        fasta_file: &str
    ) -> Self {
        let mut vargraph = VarGraph {
            multidigraph: MultiDiGraph::new(),
            nodes: HashMap::new(),
            edges: HashMap::new(),
            reference_nodes_index: BTreeMap::new(),
            variant_nodes_index: BTreeMap::new(),
            annotations: HashMap::new()
        };

        let sequence_ids: Vec<(Box<str>, usize)> = get_fasta_sequence_ids(fasta_file);
        for (sequence_id, sequence_length) in &sequence_ids {
            let reference_sequence: Box<str> = get_fasta_sequence(&*sequence_id, 1, *sequence_length, fasta_file);
            vargraph.add_reference(&*sequence_id, 1, *sequence_length, &*reference_sequence);
        }

        vargraph
    }

    pub fn from_reference_nodes(reference_nodes: Vec<&VarGraphReferenceNode>) -> Self {
        let mut vargraph = VarGraph {
            multidigraph: MultiDiGraph::new(),
            nodes: HashMap::new(),
            edges: HashMap::new(),
            reference_nodes_index: BTreeMap::new(),
            variant_nodes_index: BTreeMap::new(),
            annotations: HashMap::new()
        };

        for reference_node in reference_nodes {
            vargraph.add_reference(
                &*reference_node.chromosome,
                reference_node.start,
                reference_node.end,
                &*reference_node.sequence
            );
        }

        vargraph
    }

    pub fn add_edge_attribute(
        &mut self,
        from: usize,
        to: usize,
        key: &str,
        value: Box<dyn Any>
    ) {
        self.multidigraph.edge_exists(from, to);
        self.multidigraph.add_edge_attribute(
            from,
            to,
            key,
            value
        );
    }

    /// Adds an attribute to a node.
    ///
    /// Parameters:
    ///     id            :   Node ID.
    ///     key           :   Attribute key.
    ///     value         :   Attribute value.
    pub fn add_node_attribute(
        &mut self,
        id: usize,
        key: &str,
        value: Box<dyn Any>
    ) {
        self.multidigraph.add_node_attribute(id, key, value);
    }

    /// Add a variant.
    ///
    /// Parameters:
    /// - graph_operation   :   Graph operation.
    ///
    /// Returns:
    ///     Variant node ID.
    pub fn add_variant(
        &mut self,
        chromosome_1: &str,
        position_1: usize,
        orientation_1: VarGraphOrientations,
        strand_1: VarGraphStrands,
        chromosome_2: &str,
        position_2: usize,
        orientation_2: VarGraphOrientations,
        strand_2: VarGraphStrands,
        sequence: &str
    ) -> usize {
        let variant_node = VarGraphVariantNode::new(
            chromosome_1,
            position_1,
            orientation_1.clone(),
            strand_1.clone(),
            chromosome_2,
            position_2,
            orientation_2.clone(),
            strand_2.clone(),
            sequence
        );

        // Helper function to split a reference node at a given position
        let mut split_reference_at_position = |chromosome: &str, position: usize| -> usize {
            let ref_node_ids: Vec<usize> = self.query_reference_nodes(chromosome, position, position);
            assert_eq!(
                ref_node_ids.len(),
                1,
                "Expected to have only 1 reference node at {}:{}. Found {}.",
                chromosome,
                position,
                ref_node_ids.len()
            );
            self.split_reference_node(ref_node_ids[0], position)
        };

        // Step 1. Split the reference node at position_1
        let pos_1_ref_node_id = split_reference_at_position(&variant_node.chromosome_1, variant_node.position_1);

        // Step 2. Split the reference node at position_2
        let pos_2_ref_node_id = split_reference_at_position(&variant_node.chromosome_2, variant_node.position_2);

        // Step 3. Check if the variant node already exists
        let existing_node_ids = self.query_variant_nodes(
            &variant_node.chromosome_1,
            variant_node.position_1,
            variant_node.orientation_1.clone(),
            &variant_node.chromosome_2,
            variant_node.position_2,
            variant_node.orientation_2.clone(),
            &variant_node.sequence,
            true,
        );
        assert!(
            existing_node_ids.len() <= 1,
            "Expected to have at most 1 variant node from the query."
        );

        if let Some(&existing_node_id) = existing_node_ids.get(0) {
            // Variant node already exists
            existing_node_id
        } else {
            // Disable edges first
            let outgoing_node_ids_1: HashSet<usize>;
            let outgoing_node_ids_2: HashSet<usize>;
            if orientation_1 == VarGraphOrientations::Downstream {
                outgoing_node_ids_1 = self.get_downstream_node_ids(pos_1_ref_node_id);
            } else {
                outgoing_node_ids_1 = self.get_upstream_node_ids(pos_1_ref_node_id);
            }
            if orientation_2 == VarGraphOrientations::Downstream {
                outgoing_node_ids_2 = self.get_downstream_node_ids(pos_2_ref_node_id);
            } else {
                outgoing_node_ids_2 = self.get_upstream_node_ids(pos_2_ref_node_id);
            }
            for outgoing_node_id in outgoing_node_ids_1 {
                self.edges.get_mut(&(pos_1_ref_node_id, outgoing_node_id)).unwrap().enabled = false;
                self.edges.get_mut(&(outgoing_node_id, pos_1_ref_node_id)).unwrap().enabled = false;
            }
            for outgoing_node_id in outgoing_node_ids_2 {
                self.edges.get_mut(&(pos_2_ref_node_id, outgoing_node_id)).unwrap().enabled = false;
                self.edges.get_mut(&(outgoing_node_id, pos_2_ref_node_id)).unwrap().enabled = false;
            }

            // Add a new variant node
            let variant_node_id = self.multidigraph.add_new_node();
            self.nodes.insert(variant_node_id, variant_node.clone_box());
            let key = (
                variant_node.chromosome_1.clone(),
                variant_node.position_1,
                variant_node.orientation_1.clone(),
                variant_node.chromosome_2.clone(),
                variant_node.position_2,
                variant_node.orientation_2.clone()
            );
            self.variant_nodes_index.insert(key, variant_node_id);

            // Rewire the variation graph with the newly added variant node
            self.add_edges(
                VarGraphPort::new(pos_1_ref_node_id, strand_1.clone(), orientation_1.clone()),
                VarGraphPort::new(variant_node_id, VarGraphStrands::Forward, VarGraphOrientations::Upstream)
            );
            self.add_edges(
                VarGraphPort::new(variant_node_id, VarGraphStrands::Forward, VarGraphOrientations::Downstream),
                VarGraphPort::new(pos_2_ref_node_id, strand_2.clone(), orientation_2.clone())
            );

            variant_node_id
        }
    }

    pub fn display(&self) {
        println!("-------------------------------------------");
        println!("VarGraph ({} nodes and {} edges)", self.nodes.keys().len(), self.edges.len());
        // println!("Is acyclic: {}", self.has_eulerian_cycle());
        println!("-------------------------------------------");
        println!("Nodes");
        println!("-------------------------------------------");
        for (node_id, node) in self.nodes.iter() {
            if node.get_type() == VarGraphNodeTypes::Reference {
                println!("[REF] Node ID: {} ({}:{}-{} {})", node_id, node.get_chromosome_1(), node.get_position_1(), node.get_position_2(), node.get_sequence());
            } else {
                println!("[VAR] Node ID: {} ({})", node_id, node.get_sequence());
            }
        }
        println!("-------------------------------------------");
        println!("Edges");
        println!("-------------------------------------------");
        println!("From                 | To                  ");
        println!("-------------------------------------------");
        for ((from, to), edge) in self.edges.iter() {
            println!("{} ({}:{}) --> {} ({}:{})", from, edge.from.orientation.as_str(), edge.from.strand.as_str(), to, edge.to.orientation.as_str(), edge.to.strand.as_str());
        }
        println!("-------------------------------------------");
    }

    pub fn get_edge(&self, from: usize, to: usize) -> &VarGraphEdge {
        self.edges.get(&(from, to)).unwrap()
    }

    pub fn get_edge_attribute(&self, node_id_1: usize, node_id_2: usize, key: &str) -> Option<&Box<dyn Any>> {
        self.multidigraph.get_edge_attribute(node_id_1, node_id_2, key)
    }

    /// Get the upstream node IDs connected to 5' end of the node.
    pub fn get_upstream_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut node_ids: HashSet<usize> = HashSet::new();
        for node_id in self.multidigraph.get_outgoing_node_ids(id) {
            let edge: &VarGraphEdge = self.get_edge(id, node_id);
            if edge.from.orientation == VarGraphOrientations::Upstream {
                node_ids.insert(node_id);
            }
        }
        node_ids
    }

    /// Get the downstream node IDs connected to 3' end of the node.
    pub fn get_downstream_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut node_ids: HashSet<usize> = HashSet::new();
        for node_id in self.multidigraph.get_outgoing_node_ids(id) {
            let edge: &VarGraphEdge = self.get_edge(id, node_id);
            if edge.from.orientation == VarGraphOrientations::Downstream {
                node_ids.insert(node_id);
            }
        }
        node_ids
    }

    // /// Get the incoming node IDs from the 3' end.
    // pub fn get_incoming_3prime_node_ids(&self, id: usize) -> HashSet<usize> {
    //     let mut incoming_node_ids: HashSet<usize> = HashSet::new();
    //     for incoming_node_id in self.multidigraph.get_incoming_node_ids(id) {
    //         let edge: &VarGraphEdge = self.get_edge(id, incoming_node_id);
    //         if edge.from.orientation == VarGraphOrientations::Downstream {
    //             incoming_node_ids.insert(incoming_node_id);
    //         }
    //     }
    //     incoming_node_ids
    // }
    //
    // /// Get the incoming node IDs from the 5' end.
    // pub fn get_incoming_5prime_node_ids(&self, id: usize) -> HashSet<usize> {
    //     let mut incoming_node_ids: HashSet<usize> = HashSet::new();
    //     for incoming_node_id in self.multidigraph.get_incoming_node_ids(id) {
    //         let edge: &VarGraphEdge = self.get_edge(id, incoming_node_id);
    //         if edge.from.orientation == VarGraphOrientations::Upstream {
    //             incoming_node_ids.insert(incoming_node_id);
    //         }
    //     }
    //     incoming_node_ids
    // }
    //
    // /// Get the outgoing node IDs towards the 5' end.
    // pub fn get_outgoing_5prime_node_ids(&self, id: usize) -> HashSet<usize> {
    //     let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
    //     for outgoing_node_id in self.multidigraph.get_outgoing_node_ids(id) {
    //         let edge: &VarGraphEdge = self.get_edge(id, outgoing_node_id);
    //         if edge.from.orientation == VarGraphOrientations::Upstream {
    //             incoming_node_ids.insert(incoming_node_id);
    //         }
    //         if let Some(value) = self.multidigraph.get_edge_attribute(id, outgoing_node_id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
    //             if let Some(direction) = value.downcast_ref::<Box<str>>() {
    //                 if &**direction == VarGraphEdgeDirections::ThreeToFivePrime.as_str() {
    //                     outgoing_node_ids.insert(outgoing_node_id);
    //                 }
    //             }
    //         }
    //     }
    //     outgoing_node_ids
    // }
    //
    // /// Get the outgoing node IDs towards the 3' end.
    // pub fn get_outgoing_3prime_node_ids(&self, id: usize) -> HashSet<usize> {
    //     let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
    //     for outgoing_node_id in self.multidigraph.get_outgoing_node_ids(id) {
    //         if let Some(value) = self.multidigraph.get_edge_attribute(id, outgoing_node_id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
    //             if let Some(direction) = value.downcast_ref::<Box<str>>() {
    //                 if &**direction == VarGraphEdgeDirections::FiveToThreePrime.as_str() {
    //                     outgoing_node_ids.insert(outgoing_node_id);
    //                 }
    //             }
    //         }
    //     }
    //     outgoing_node_ids
    // }

    /// Fetch a VarGraphNode object.
    ///
    /// Parameters:
    ///     id        :   Node ID.
    ///
    /// Returns:
    /// - Reference to a VarGraphNode object.
    pub fn get_node(&self, id: usize) -> &dyn VarGraphNode {
        match self.nodes.get(&id) {
            Some(value) => value.as_ref(),
            None => {
                panic!("A node with ID {} does not exist.", id);
            }
        }
    }

    /// Get the total number of nodes.
    ///
    /// Returns:
    /// - Total number of nodes.
    pub fn get_nodes_count(&self) -> usize {
        self.nodes.keys().len()
    }

    pub fn get_node_type(&self, id: usize) -> VarGraphNodeTypes {
        self.nodes.get(&id).unwrap().get_type()
    }

    /// Get reference node Ids.
    pub fn get_reference_node_ids(&self) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        for (key, value) in self.reference_nodes_index.iter() {
            node_ids.push(*value);
        }
        node_ids
    }

    /// Get variant node Ids.
    pub fn get_variant_node_ids(&self) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        for (key, value) in self.variant_nodes_index.iter() {
            node_ids.push(*value);
        }
        node_ids
    }

    pub fn get_linearized_contigs(&self, variant_node_ids: Vec<usize>) -> Vec<VarGraphPath> {
        #[derive(Debug,Clone)]
        struct TraversalState {
            visited: HashSet<usize>,
            path: VecDeque<usize>
        }

        // Step 1. Identify node paths
        let mut node_paths: HashSet<VecDeque<usize>> = HashSet::new();
        for variant_node_id in variant_node_ids.iter() {
            let mut states: Vec<TraversalState> = Vec::new();
            let mut queue: VecDeque<TraversalState> = VecDeque::new();

            // Create a new path from the variant node
            let initial_state: TraversalState = TraversalState {
                visited: HashSet::from_iter(vec![*variant_node_id]),
                path: VecDeque::from_iter(vec![*variant_node_id])
            };

            // Traverse in the downstream direction (5' to 3') of the variant node
            for node_id in self.get_downstream_node_ids(*variant_node_id).iter() {
                let mut state: TraversalState = initial_state.clone();
                state.visited.insert(*node_id);
                state.path.push_back(*node_id);
                queue.push_back(state);
            }
            while let Some(mut state) = queue.pop_front() {
                let last_node_id: usize = *state.path.back().unwrap();

                // Get next node IDs that have enabled edges from the last node
                let mut next_node_ids: HashSet<usize> = HashSet::new();
                for outgoing_node_id in self.get_upstream_node_ids(last_node_id) {
                    if self.get_edge(last_node_id, outgoing_node_id).enabled && state.visited.contains(&outgoing_node_id) == false {
                        next_node_ids.insert(outgoing_node_id);
                    }
                }
                for outgoing_node_id in self.get_downstream_node_ids(last_node_id) {
                    if self.get_edge(last_node_id, outgoing_node_id).enabled && state.visited.contains(&outgoing_node_id) == false {
                        next_node_ids.insert(outgoing_node_id);
                    }
                }

                if next_node_ids.is_empty() {
                    states.push(state.clone());
                } else {
                    for next_node_id in next_node_ids {
                        let mut new_state: TraversalState = state.clone();
                        new_state.path.push_back(next_node_id);
                        new_state.visited.insert(next_node_id);
                        queue.push_back(new_state);
                    }
                }
            }

            // Traverse in the upstream direction (3' to 5') of the variant node
            queue.clear();
            for state_ in states.iter() {
                for node_id in self.get_upstream_node_ids(*variant_node_id).iter() {
                    let mut state: TraversalState = state_.clone();
                    state.path.push_front(*node_id);
                    state.visited.insert(*node_id);
                    queue.push_back(state);
                }
            }
            states.clear();
            while let Some(mut state) = queue.pop_front() {
                let first_node_id: usize = *state.path.front().unwrap();

                // Get next node IDs that have enabled edges from the first node
                let mut next_node_ids: HashSet<usize> = HashSet::new();
                for outgoing_node_id in self.get_upstream_node_ids(first_node_id) {
                    if self.get_edge(first_node_id, outgoing_node_id).enabled && state.visited.contains(&outgoing_node_id) == false {
                        next_node_ids.insert(outgoing_node_id);
                    }
                }
                for outgoing_node_id in self.get_downstream_node_ids(first_node_id) {
                    if self.get_edge(first_node_id, outgoing_node_id).enabled && state.visited.contains(&outgoing_node_id) == false {
                        next_node_ids.insert(outgoing_node_id);
                    }
                }

                if next_node_ids.is_empty() {
                    states.push(state.clone());
                } else {
                    for next_node_id in next_node_ids {
                        let mut new_state: TraversalState = state.clone();
                        new_state.path.push_front(next_node_id);
                        new_state.visited.insert(next_node_id);
                        queue.push_back(new_state);
                    }
                }
            }

            for state in states {
                node_paths.insert(state.path);
            }
        }

        // Step 2. Convert node paths to segments
        let mut paths: Vec<VarGraphPath> = Vec::new();
        for node_path in node_paths.iter() {
            let mut path: VarGraphPath = VarGraphPath::new();

            // First node
            let first_node_id: usize = *node_path.get(0).unwrap();
            let second_node_id: usize = *node_path.get(1).unwrap();
            let first_node: &dyn VarGraphNode = self.get_node(first_node_id);
            let edge: &VarGraphEdge = self.get_edge(first_node_id, second_node_id);
            let end: usize = (first_node.get_sequence_length() as usize).saturating_sub(1);
            let segment = VarGraphSegment::new(
                first_node_id,
                0,
                end,
                &*first_node.get_sequence(),
                None,
                Some(edge.from.clone())
            );
            path.push_back(segment);

            // Second node to the second last node
            for i in 1..=node_path.len() - 2 {
                let prev_node_id: usize = *node_path.get(i - 1).unwrap();
                let curr_node_id: usize = *node_path.get(i).unwrap();
                let next_node_id: usize = *node_path.get(i + 1).unwrap();
                let curr_node: &dyn VarGraphNode = self.get_node(curr_node_id);
                let prev_edge: &VarGraphEdge = self.get_edge(prev_node_id, curr_node_id);
                let next_edge: &VarGraphEdge = self.get_edge(curr_node_id, next_node_id);
                let end: usize = (curr_node.get_sequence_length() as usize).saturating_sub(1);
                // let sequence: Box<str> = if prev_edge.to.strand == VarGraphStrands::Forward {
                //     curr_node.get_sequence().clone()
                // } else {
                //     reverse_complement(&*curr_node.get_sequence())
                // };
                let segment = VarGraphSegment::new(
                    curr_node_id,
                    0,
                    end,
                    &*curr_node.get_sequence(),
                    Some(prev_edge.to.clone()),
                    Some(next_edge.from.clone())
                );
                path.push_back(segment);
            }

            // Last node
            let last_node_id: usize = *node_path.get(node_path.len() - 1).unwrap();
            let second_last_node_id: usize = *node_path.get(node_path.len() - 2).unwrap();
            let last_node: &dyn VarGraphNode = self.get_node(last_node_id);
            let edge: &VarGraphEdge = self.get_edge(second_last_node_id, last_node_id);
            let end: usize = (last_node.get_sequence_length() as usize).saturating_sub(1);
            let segment = VarGraphSegment::new(
                last_node_id,
                0,
                end,
                &*last_node.get_sequence(),
                Some(edge.to.clone()),
                None
            );
            path.push_back(segment);

            paths.push(path);
        }

        paths
    }

    // pub fn get_all_linearized_contigs(&self) -> Vec<VarGraphPath> {
    //     // Step 1. Get all linearized contigs involving all variant nodes
    //     let mut vargraph_paths: Vec<VarGraphPath> = self.get_linearized_contigs(self.get_variant_node_ids());
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

    /// Find reference node IDs for a query region.
    ///
    /// Parameters:
    /// - chromosome    :   Query chromosome.
    /// - start         :   Query start.
    /// - end           :   Query end.
    pub fn query_reference_nodes(
        &self,
        chromosome: &str,
        start: usize,
        end: usize
    ) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        let lower_bound: (Box<str>, usize, usize) = (
            chromosome.into(),
            0,
            0
        );
        let upper_bound: (Box<str>, usize, usize) = (
            chromosome.into(),
            usize::MAX,
            usize::MAX
        );
        for ((cmp_chromosome,cmp_position_1,cmp_position_2), value)
            in self.reference_nodes_index.range((Included(lower_bound), Included(upper_bound))) {
            if &**cmp_chromosome == chromosome &&
                *cmp_position_1 <= start &&
                *cmp_position_2 >= end {
                node_ids.push(*value);
            }
        }
        node_ids
    }

    /// Find variant node IDs for a query region.
    ///
    /// Parameters:
    /// - chromosome_1      :   Query chromosome 1.
    /// - position_1        :   Query position 1.
    /// - orientation_1     :   Query orientation 1.
    /// - chromosome_2      :   Query chromosome 2.
    /// - position_2        :   Query position 2.
    /// - orientation_2     :   Query orientation 2.
    /// - sequence          :   Query sequence.
    /// - match_sequence    :   If true, the sequences must match.
    pub fn query_variant_nodes(
        &self,
        chromosome_1: &str,
        position_1: usize,
        orientation_1: VarGraphOrientations,
        chromosome_2: &str,
        position_2: usize,
        orientation_2: VarGraphOrientations,
        sequence: &str,
        match_sequence: bool
    ) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        let lower_bound: (Box<str>,usize,VarGraphOrientations,Box<str>,usize,VarGraphOrientations) = (
            chromosome_1.into(),
            0,
            orientation_1.clone(),
            chromosome_2.into(),
            0,
            orientation_2.clone()
        );
        let upper_bound: (Box<str>,usize,VarGraphOrientations,Box<str>,usize,VarGraphOrientations) = (
            chromosome_1.into(),
            usize::MAX,
            orientation_1.clone(),
            chromosome_2.into(),
            usize::MAX,
            orientation_2.clone()
        );
        for ((cmp_chromosome_1,cmp_position_1,cmp_orientation_1,cmp_chromosome_2,cmp_position_2,cmp_orientation_2), value)
            in self.variant_nodes_index.range((Included(lower_bound), Included(upper_bound))) {
            if &**cmp_chromosome_1 != chromosome_1 ||
                &**cmp_chromosome_2 != chromosome_2 ||
                *cmp_position_1 > position_1 ||
                *cmp_position_2 < position_2 ||
                *cmp_orientation_1 != orientation_1 ||
                *cmp_orientation_2 != orientation_2 {
                continue;
            }
            let cmp_sequence = self.get_node(*value).get_sequence();
            let sequence_match = !match_sequence || &*cmp_sequence == sequence;
            if sequence_match {
                node_ids.push(*value);
            }
        }
        node_ids
    }
}

// Private methods
impl VarGraph {

    /// Add edges (5' to 3' and 3' to 5') between two nodes.
    ///
    /// Parameters:
    /// - from      :   VarGraph Port object.
    /// - to        :   VarGraph Port object.
    fn add_edges(
        &mut self,
        from: VarGraphPort,
        to: VarGraphPort,
    ) {
        self.multidigraph.add_edge(from.node_id, to.node_id);
        self.multidigraph.add_edge(to.node_id, from.node_id);
        self.edges.insert((from.node_id, to.node_id), VarGraphEdge::new(from.clone(), to.clone()));
        self.edges.insert((to.node_id, from.node_id), VarGraphEdge::new(to, from));
    }

    /// Adds a reference.
    ///
    /// Parameters:
    /// - chromosome      :   Chromosome.
    /// - start           :   Start.
    /// - end             :   End.
    /// - sequence        :   Sequence (forward/coding strand sequence).
    ///
    /// Returns:
    /// - Reference node ID.
    fn add_reference(
        &mut self,
        chromosome: &str,
        start: usize,
        end: usize,
        sequence: &str
    ) -> usize {
        // Step 1. Error checking
        if start > end {
            panic!("'start' must be equal to or smaller than 'end'.");
        }
        if sequence.len() != (end - start + 1) {
            panic!("'sequence' length must match the length from 'start' to 'end'.");
        }

        // Step 2. Add the new reference node
        let new_node_id: usize = self.multidigraph.add_new_node();
        let reference_node: VarGraphReferenceNode = VarGraphReferenceNode::new(
            chromosome, start, end, sequence
        );
        self.nodes.insert(new_node_id, Box::new(reference_node));

        // Step 4. Index the new reference node
        let key: (Box<str>, usize, usize) = (
            chromosome.to_string().into_boxed_str(),
            start,
            end
        );
        self.reference_nodes_index.insert(key, new_node_id);

        new_node_id
    }

    /// Removes a node.
    ///
    /// Parameters:
    /// - id    :   Node ID.
    fn remove_node(&mut self, id: usize) {
        if self.nodes.get(&id).unwrap().get_type() == VarGraphNodeTypes::Reference {
            let reference_node: &VarGraphReferenceNode = self.get_node(id).as_any().downcast_ref::<VarGraphReferenceNode>().unwrap();
            let index_key: (Box<str>, usize, usize) = (
                reference_node.chromosome.clone(),
                reference_node.start,
                reference_node.end
            );
            self.reference_nodes_index.remove(&index_key);
        } else {
            let variant_node: &VarGraphVariantNode = self.get_node(id).as_any().downcast_ref::<VarGraphVariantNode>().unwrap();
            let index_key: (Box<str>,usize,VarGraphOrientations,Box<str>,usize,VarGraphOrientations) = (
                variant_node.chromosome_1.clone(),
                variant_node.position_1,
                variant_node.orientation_1.clone(),
                variant_node.chromosome_2.clone(),
                variant_node.position_2,
                variant_node.orientation_2.clone()
            );
            self.variant_nodes_index.remove(&index_key);
        }
        self.nodes.remove(&id);
        self.edges.retain(|key, _| key.0 != id && key.1 != id);
        self.multidigraph.remove_node(id);
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
    fn split_reference_node(&mut self, id: usize, position: usize) -> usize {
        assert_eq!(self.nodes.contains_key(&id), true, "A node with ID {} does not exist.", id);

        // Step 1. Fetch the reference node object
        let node: &dyn VarGraphNode = self.get_node(id);
        let reference_node: &VarGraphReferenceNode = node.as_any().downcast_ref::<VarGraphReferenceNode>().unwrap();
        if reference_node.sequence.len() <= 1 {
            // The node is already independent
            return id;
        }

        // Step 2. Get all outgoing node IDs to 5' and 3' ends from the reference node
        let upstream_node_ids: HashSet<usize> = self.get_upstream_node_ids(id);
        let downstream_node_ids: HashSet<usize> = self.get_downstream_node_ids(id);

        // Step 3. Split the reference node
        let returning_node_id: usize;
        if position == reference_node.start {
            // Create 2 new VarGraphReferenceNode objects
            // At position 'position'
            let node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                position,
                position,
                &reference_node.sequence[0..1]
            );

            // From (position+1) to the end of the reference node
            let node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                position + 1,
                reference_node.end,
                &reference_node.sequence[1..]
            );

            // Generate new nodes
            let node_id_1 = self.multidigraph.add_new_node();
            let node_id_2 = self.multidigraph.add_new_node();

            // Add the new nodes
            self.nodes.insert(node_id_1.clone(), Box::new(node_1.clone()));
            self.nodes.insert(node_id_2.clone(), Box::new(node_2.clone()));

            // Add the new nodes to the index
            let key_1: (Box<str>, usize, usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1);
            let key_2: (Box<str>, usize, usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in upstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_1, existing_edge.from.strand.clone(), VarGraphOrientations::Upstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Downstream),
                );
            }

            // Add an edge between the new nodes
            self.add_edges(
                VarGraphPort::new(node_id_1, VarGraphStrands::Forward, VarGraphOrientations::Downstream),
                VarGraphPort::new(node_id_2, VarGraphStrands::Forward, VarGraphOrientations::Upstream),
            );

            // Add edges from all original outgoing nodes towards the 3' end from node_2
            for node_id in downstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_2, existing_edge.from.strand.clone(), VarGraphOrientations::Downstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Upstream),
                );
            }

            // Remove the old node
            self.remove_node(id);

            // Assign the returning node ID
            returning_node_id = node_id_1;
        } else if position == reference_node.end {
            // Create 2 new VarGraphReferenceNode objects
            // From start of the reference node to (position-1)
            let node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                reference_node.start,
                position - 1,
                &reference_node.sequence[0..reference_node.sequence.len() - 1]
            );

            // At position 'position'
            let node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                position,
                position,
                &reference_node.sequence[reference_node.sequence.len() - 1..]
            );

            // Generate new nodes
            let node_id_1 = self.multidigraph.add_new_node();
            let node_id_2 = self.multidigraph.add_new_node();

            // Add the new nodes
            self.nodes.insert(node_id_1.clone(), Box::new(node_1.clone()));
            self.nodes.insert(node_id_2.clone(), Box::new(node_2.clone()));

            // Add the new nodes to the index
            let key_1: (Box<str>, usize, usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1);
            let key_2: (Box<str>, usize, usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in upstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_1, existing_edge.from.strand.clone(), VarGraphOrientations::Upstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Downstream)
                );
            }

            // Add an edge between the new nodes
            self.add_edges(
                VarGraphPort::new(node_id_1, VarGraphStrands::Forward, VarGraphOrientations::Downstream),
                VarGraphPort::new(node_id_2, VarGraphStrands::Forward, VarGraphOrientations::Upstream)
            );

            // Add edges from all original outgoing nodes towards the 3' end from node_2
            for node_id in downstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_2, existing_edge.from.strand.clone(), VarGraphOrientations::Downstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Upstream)
                );
            }

            // Remove the old node
            self.remove_node(id);

            // Assign the returning node ID
            returning_node_id = node_id_2;
        } else {
            // Create 3 new VarGraphReferenceNode objects
            // From the start of the reference node to (position-1)
            let node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                reference_node.start,
                position - 1,
                &reference_node.sequence[..(position - reference_node.start)]
            );

            // At position 'position'
            let node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                position,
                position,
                &reference_node.sequence[(position - reference_node.start)..(position - reference_node.start + 1)]
            );

            // From (position+1) to the end of the reference node
            let node_3: VarGraphReferenceNode = VarGraphReferenceNode::new(
                &*reference_node.chromosome,
                position + 1,
                reference_node.end,
                &reference_node.sequence[(position - reference_node.start + 1)..]
            );

            // Generate new nodes
            let node_id_1 = self.multidigraph.add_new_node();
            let node_id_2 = self.multidigraph.add_new_node();
            let node_id_3 = self.multidigraph.add_new_node();

            // Add the new nodes
            self.nodes.insert(node_id_1.clone(), Box::new(node_1.clone()));
            self.nodes.insert(node_id_2.clone(), Box::new(node_2.clone()));
            self.nodes.insert(node_id_3.clone(), Box::new(node_3.clone()));

            // Add the new nodes to the index
            let key_1: (Box<str>, usize, usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1.clone());
            let key_2: (Box<str>, usize, usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2.clone());
            let key_3: (Box<str>, usize, usize) = (
                node_3.chromosome.clone(),
                node_3.start,
                node_3.end
            );
            self.reference_nodes_index.insert(key_3, node_id_3.clone());

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in upstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_1, existing_edge.from.strand.clone(), VarGraphOrientations::Upstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Downstream)
                );
            }

            // Add an edge between the new nodes
            self.add_edges(
                VarGraphPort::new(node_id_1, VarGraphStrands::Forward, VarGraphOrientations::Downstream),
                VarGraphPort::new(node_id_2, VarGraphStrands::Forward, VarGraphOrientations::Upstream)
            );
            self.add_edges(
                VarGraphPort::new(node_id_2, VarGraphStrands::Forward, VarGraphOrientations::Downstream),
                VarGraphPort::new(node_id_3, VarGraphStrands::Forward, VarGraphOrientations::Upstream)
            );

            // Add edges from all original outgoing nodes towards the 3' end from node_3
            for node_id in downstream_node_ids.iter() {
                let existing_edge: &VarGraphEdge = self.edges.get(&(id, *node_id)).unwrap();
                self.add_edges(
                    VarGraphPort::new(node_id_3, existing_edge.from.strand.clone(), VarGraphOrientations::Downstream),
                    VarGraphPort::new(*node_id, existing_edge.to.strand.clone(), VarGraphOrientations::Upstream)
                );
            }

            // Remove the old node
            self.remove_node(id);

            // Assign the returning node ID
            returning_node_id = node_id_2;
        }

        returning_node_id
    }
}
