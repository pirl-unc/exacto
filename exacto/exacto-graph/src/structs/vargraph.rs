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


//! **vargraph** implements a variation graph (a directed graph of sequences that can store multiple edges).
//!
//! # Examples
//!
//! ## Construct a VarGraph (a directed graph of sequences that can store multiple edges)
//!
//! ```no_run
//! use exacto_graph::structs::vargraph::*;
//!
//! /// Variation graph
//! let mut vargraph = VarGraph::new();
//! vargraph.add_reference("chr1", 1, 10, "ACGTACGATG");
//! ```


use std::any::Any;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::Bound::Included;

use crate::common::constants::*;
use crate::common::parsers::parse_graph_operation;
use crate::structs::multidigraph::MultiDiGraph;
use crate::structs::vargraph_variant_node::VarGraphVariantNode;
use crate::structs::vargraph_reference_node::VarGraphReferenceNode;
use crate::structs::vargraph_segment::VarGraphSegment;
use crate::structs::vargraph_path::VarGraphPath;
use crate::traits::vargraph_node::VarGraphNode;


#[derive(Debug)]
pub struct VarGraph {
    pub multidigraph: MultiDiGraph,

    /// Nodes in the variation graph.
    ///
    /// HashMap
    ///     Key   :   Node ID.
    ///     Value :   SequenceGraphNode object.
    pub nodes: HashMap<usize,Box<dyn VarGraphNode>>,

    /// Reference nodes index.
    ///
    /// BTreeMap:
    ///     Key     :   (chromosome_1, position_1, chromosome_2, position_2).
    ///     Value   :   Node ID.
    pub reference_nodes_index: BTreeMap<(Box<str>,usize,Box<str>,usize),usize>,

    /// Variant nodes index.
    ///
    /// BTreeMap:
    ///     Key     :   (chromosome_1, position_1, orientation_1, chromosome_2, position_2, orientation_2).
    ///     Value   :   Node ID.
    pub variant_nodes_index: BTreeMap<(Box<str>,usize,VarGraphOrientations,Box<str>,usize,VarGraphOrientations),usize>,

    /// Region annotations.
    ///
    /// Outer HashMap
    ///     Key     : (chromosome, start, end).
    ///     Value   : An inner HashMap representing the region gene_annotation.
    ///
    /// Inner HashMap
    ///     Key     : Attribute name.
    ///     Value   : Attribute value.
    pub annotations: HashMap<(Box<str>,usize,usize),HashMap<Box<str>,Box<dyn Any>>>,
}

impl VarGraph {
    /// Constructor
    pub fn new() -> Self {
        VarGraph {
            multidigraph: MultiDiGraph::new(),
            nodes: HashMap::new(),
            reference_nodes_index: BTreeMap::new(),
            variant_nodes_index: BTreeMap::new(),
            annotations: HashMap::new()
        }
    }

    /// Add edges (5' to 3' and 3' to 5') between two nodes.
    ///
    /// Parameters:
    /// - node_id_1     :   Node ID 1.
    /// - node_id_2     :   Node ID 2.
    /// - direction     :   Direction (node_id_1 to node_id_2).
    fn add_edges(
        &mut self,
        node_id_1: usize,
        node_id_2: usize,
        direction: VarGraphEdgeDirections
    ) {
        self.multidigraph.add_edge(node_id_1, node_id_2);
        self.multidigraph.add_edge(node_id_2, node_id_1);
        let direction_1: Box<str>;
        let direction_2: Box<str>;
        if direction == VarGraphEdgeDirections::FiveToThreePrime {
            direction_1 = VarGraphEdgeDirections::FiveToThreePrime.as_str().into();
            direction_2 = VarGraphEdgeDirections::ThreeToFivePrime.as_str().into();
        } else {
            direction_1 = VarGraphEdgeDirections::ThreeToFivePrime.as_str().into();
            direction_2 = VarGraphEdgeDirections::FiveToThreePrime.as_str().into();
        }
        self.multidigraph.add_edge_attribute(
            node_id_1,
            node_id_2,
            VarGraphEdgeAttributeKeys::Direction.as_str(),
            Box::new(direction_1)
        );
        self.multidigraph.add_edge_attribute(
            node_id_2,
            node_id_1,
            VarGraphEdgeAttributeKeys::Direction.as_str(),
            Box::new(direction_2)
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
    pub fn add_reference(
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
        let key: (Box<str>,usize,Box<str>,usize) = (
            chromosome.to_string().into_boxed_str(),
            start,
            chromosome.to_string().into_boxed_str(),
            end
        );
        self.reference_nodes_index.insert(key, new_node_id);

        new_node_id
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
        graph_operation: &str
    ) -> usize {
        let variant_node: VarGraphVariantNode = parse_graph_operation(graph_operation);

        // Helper function to split a reference node at a given position
        let mut split_reference_at_position = |chromosome: &str, position: usize| -> usize {
            let ref_node_ids: Vec<usize> = self.query_reference_nodes(chromosome, position, chromosome, position);
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

            // Helper function to add edges
            let mut add_variant_edges = |ref_node_id: usize, orientation: VarGraphOrientations| {
                let direction = if orientation == VarGraphOrientations::Upstream {
                    VarGraphEdgeDirections::ThreeToFivePrime
                } else {
                    VarGraphEdgeDirections::FiveToThreePrime
                };
                self.add_edges(ref_node_id, variant_node_id, direction);
            };

            // Rewire the variation graph with the newly added variant node
            add_variant_edges(pos_1_ref_node_id, variant_node.orientation_1.clone());
            add_variant_edges(pos_2_ref_node_id, variant_node.orientation_2.clone());

            variant_node_id
        }
    }

    /// Get node IDs incoming from the 3' end.
    pub fn get_incoming_3prime_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut incoming_node_ids: HashSet<usize> = HashSet::new();
        for incoming_node_id in self.multidigraph.get_incoming_node_ids(id) {
            if let Some(value) = self.multidigraph.get_edge_attribute(incoming_node_id, id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
                if let Some(direction) = value.downcast_ref::<Box<str>>() {
                    if &**direction == VarGraphEdgeDirections::ThreeToFivePrime.as_str() {
                        incoming_node_ids.insert(incoming_node_id);
                    }
                }
            }
        }
        incoming_node_ids
    }

    /// Get node IDs incoming from the 5' end.
    pub fn get_incoming_5prime_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut incoming_node_ids: HashSet<usize> = HashSet::new();
        for incoming_node_id in self.multidigraph.get_incoming_node_ids(id) {
            if let Some(value) = self.multidigraph.get_edge_attribute(incoming_node_id, id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
                if let Some(direction) = value.downcast_ref::<Box<str>>() {
                    if &**direction == VarGraphEdgeDirections::FiveToThreePrime.as_str() {
                        incoming_node_ids.insert(incoming_node_id);
                    }
                }
            }
        }
        incoming_node_ids
    }

    /// Get node IDs outgoing towards the 5' end.
    pub fn get_outgoing_5prime_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
        for outgoing_node_id in self.multidigraph.get_outgoing_node_ids(id) {
            if let Some(value) = self.multidigraph.get_edge_attribute(id, outgoing_node_id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
                if let Some(direction) = value.downcast_ref::<Box<str>>() {
                    if &**direction == VarGraphEdgeDirections::ThreeToFivePrime.as_str() {
                        outgoing_node_ids.insert(outgoing_node_id);
                    }
                }
            }
        }
        outgoing_node_ids
    }

    /// Get node IDs outgoing towards the 3' end.
    pub fn get_outgoing_3prime_node_ids(&self, id: usize) -> HashSet<usize> {
        let mut outgoing_node_ids: HashSet<usize> = HashSet::new();
        for outgoing_node_id in self.multidigraph.get_outgoing_node_ids(id) {
            if let Some(value) = self.multidigraph.get_edge_attribute(id, outgoing_node_id, VarGraphEdgeAttributeKeys::Direction.as_str()) {
                if let Some(direction) = value.downcast_ref::<Box<str>>() {
                    if &**direction == VarGraphEdgeDirections::FiveToThreePrime.as_str() {
                        outgoing_node_ids.insert(outgoing_node_id);
                    }
                }
            }
        }
        outgoing_node_ids
    }

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

    /// Find reference node IDs for a query region.
    ///
    /// Parameters:
    /// - chromosome_1      :   Query chromosome 1.
    /// - position_1        :   Query position 1.
    /// - chromosome_2      :   Query chromosome 2.
    /// - position_2        :   Query position 2.
    pub fn query_reference_nodes(
        &self,
        chromosome_1: &str,
        position_1: usize,
        chromosome_2: &str,
        position_2: usize,
    ) -> Vec<usize> {
        let mut node_ids: Vec<usize> = Vec::new();
        let lower_bound: (Box<str>,usize,Box<str>,usize) = (
            chromosome_1.into(),
            0,
            chromosome_2.into(),
            0
        );
        let upper_bound: (Box<str>,usize,Box<str>,usize) = (
            chromosome_1.into(),
            usize::MAX,
            chromosome_2.into(),
            usize::MAX
        );
        for ((cmp_chromosome_1,cmp_position_1,cmp_chromosome_2,cmp_position_2), value)
            in self.reference_nodes_index.range((Included(lower_bound), Included(upper_bound))) {
            if &**cmp_chromosome_1 == chromosome_1 &&
                &**cmp_chromosome_2 == chromosome_2 &&
                *cmp_position_1 <= position_1 &&
                *cmp_position_2 >= position_2 {
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

    /// Removes a node.
    ///
    /// Parameters:
    /// - id    :   Node ID.
    pub fn remove_node(&mut self, id: usize) {
        if self.nodes.get(&id).unwrap().get_type() == VarGraphNodeTypes::Reference {
            let reference_node: &VarGraphReferenceNode = self.get_node(id).as_any().downcast_ref::<VarGraphReferenceNode>().unwrap();
            let index_key: (Box<str>,usize,Box<str>,usize) = (
                reference_node.chromosome.clone(),
                reference_node.start,
                reference_node.chromosome.clone(),
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
        self.multidigraph.remove_node(id);
    }

    /// Splits a reference node. The result is that there will be a reference node at
    /// position 'position' with sequence length of 1.
    ///
    /// Parameters:
    /// - node_id   :   Reference node ID.
    /// - position  :   Position to split.
    ///
    /// Returns:
    /// Node ID at position 'position'.
    pub fn split_reference_node(&mut self, id: usize, position: usize) -> usize {
        assert_eq!(self.nodes.contains_key(&id), true, "A node with ID {} does not exist.", id);

        // Step 1. Fetch the reference node object
        let node: &dyn VarGraphNode = self.get_node(id);
        let reference_node: &VarGraphReferenceNode = node.as_any().downcast_ref::<VarGraphReferenceNode>().unwrap();
        if reference_node.sequence.len() <= 1 {
            // The node is already independent
            return id;
        }

        // Step 2. Get all outgoing node IDs to 5' and 3' ends from the reference node
        let outgoing_5prime_node_ids: HashSet<usize> = self.get_outgoing_5prime_node_ids(id);
        let outgoing_3prime_node_ids: HashSet<usize> = self.get_outgoing_3prime_node_ids(id);

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
            let key_1: (Box<str>,usize,Box<str>,usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.chromosome.clone(),
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1);
            let key_2: (Box<str>,usize,Box<str>,usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.chromosome.clone(),
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in outgoing_5prime_node_ids.iter() {
                self.add_edges(node_id_1, *node_id, VarGraphEdgeDirections::ThreeToFivePrime);
            }

            // Add an edge between the new nodes
            self.add_edges(node_id_1, node_id_2, VarGraphEdgeDirections::FiveToThreePrime);

            // Add edges from all original outgoing nodes towards the 3' end from node_2
            for node_id in outgoing_3prime_node_ids.iter() {
                self.add_edges(node_id_2, *node_id, VarGraphEdgeDirections::FiveToThreePrime);
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
            let key_1: (Box<str>,usize,Box<str>,usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.chromosome.clone(),
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1);
            let key_2: (Box<str>,usize,Box<str>,usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.chromosome.clone(),
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2);

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in outgoing_5prime_node_ids.iter() {
                self.add_edges(node_id_1, *node_id, VarGraphEdgeDirections::ThreeToFivePrime);
            }

            // Add an edge between the new nodes
            self.add_edges(node_id_1, node_id_2, VarGraphEdgeDirections::FiveToThreePrime);

            // Add edges from all original outgoing nodes towards the 3' end from node_2
            for node_id in outgoing_3prime_node_ids.iter() {
                self.add_edges(node_id_2, *node_id, VarGraphEdgeDirections::FiveToThreePrime);
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
            let key_1: (Box<str>,usize,Box<str>,usize) = (
                node_1.chromosome.clone(),
                node_1.start,
                node_1.chromosome.clone(),
                node_1.end
            );
            self.reference_nodes_index.insert(key_1, node_id_1.clone());
            let key_2: (Box<str>,usize,Box<str>,usize) = (
                node_2.chromosome.clone(),
                node_2.start,
                node_2.chromosome.clone(),
                node_2.end
            );
            self.reference_nodes_index.insert(key_2, node_id_2.clone());
            let key_3: (Box<str>,usize,Box<str>,usize) = (
                node_3.chromosome.clone(),
                node_3.start,
                node_3.chromosome.clone(),
                node_3.end
            );
            self.reference_nodes_index.insert(key_3, node_id_3.clone());

            // Add edges from all original outgoing nodes towards the 5' end from node_1
            for node_id in outgoing_5prime_node_ids.iter() {
                self.add_edges(node_id_1, *node_id, VarGraphEdgeDirections::ThreeToFivePrime);
            }

            // Add an edge between the new nodes
            self.add_edges(node_id_1, node_id_2, VarGraphEdgeDirections::FiveToThreePrime);
            self.add_edges(node_id_2, node_id_3, VarGraphEdgeDirections::FiveToThreePrime);

            // Add edges from all original outgoing nodes towards the 3' end from node_3
            for node_id in outgoing_3prime_node_ids.iter() {
                self.add_edges(node_id_3, *node_id, VarGraphEdgeDirections::FiveToThreePrime);
            }

            // Remove the old node
            self.remove_node(id);

            // Assign the returning node ID
            returning_node_id = node_id_2;
        }

        returning_node_id
    }
}