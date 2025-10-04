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


//! **graph** implements an undirected multi graph.
//!
//! The **graph** struct serves as a base abstraction.
//!
//! # Examples
//!
//! ## Construct a Graph (undirected multi graph)
//!
//! ```no_run
//! # use std::io;
//! use exacto_graph::structs::graph::*;
//!
//! /// Undirected multi graph
//! let mut graph = Graph::new();
//! let node_1: usize = graph.add_new_node();
//! let node_2: usize = graph.add_new_node();
//! let node_3: usize = graph.add_new_node();
//! let node_4: usize = graph.add_new_node();
//! graph.add_edge(node_1, node_2);
//! graph.add_edge(node_2, node_3);
//! graph.add_edge(node_3, node_4);
//! graph.add_edge(node_1, node_3);
//! graph.remove_node(node_2);
//!
//! /// The in-degree for the node_3 is 2 (1-->, D-->C)
//! let result: usize = graph.get_in_degree(node_3);
//! ```


use exacto_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::any::Any;


/// Graph
#[derive(Debug)]
pub struct Graph {
    /// Nodes in the graph.
    pub nodes: HashSet<usize>,

    /// Outgoing edges of the graph.
    ///
    /// - HashMap
    ///     Key     :   Source node ID.
    ///     Value   :   HashSet of destination node IDs.
    pub outgoing_edges: HashMap<usize, HashSet<usize>>,

    /// Incoming edges of the graph.
    ///
    /// HashMap
    ///     Key     :   Destination node ID.
    ///     Value   :   HashSet of source node IDs.
    pub incoming_edges: HashMap<usize, HashSet<usize>>,

    /// Node attributes.
    ///
    /// Outer HashMap
    ///     Key     :   Node ID.
    ///     Value   :   HashMap of node attributes.
    ///
    /// Inner HashMap
    ///     Key     :   Attribute name.
    ///     Value   :   Attribute value.
    pub node_attributes: HashMap<usize,HashMap<Box<str>,Box<dyn Any>>>,

    /// Edge attributes.
    ///
    /// Outer HashMap
    ///     Key     :   Tuple (source node ID, destination node ID).
    ///     Value   :   HashMap representing the node attributes.
    ///
    /// Inner HashMap
    ///     Key     :   Attribute name.
    ///     Value   :   Attribute value.
    pub edge_attributes: HashMap<(usize,usize),HashMap<Box<str>,Box<dyn Any>>>,

    node_id_counter: usize
}

impl Graph {

    pub fn new() -> Self {
        Graph {
            nodes: HashSet::new(),
            outgoing_edges: HashMap::new(),
            incoming_edges: HashMap::new(),
            node_attributes: HashMap::new(),
            edge_attributes: HashMap::new(),
            node_id_counter: 1
        }
    }

    /// Add an edge in both directions for an undirected graph.
    /// The nodes 'from' and 'to' must be in the graph already.
    ///
    /// Parameters:
    /// - from    :   Node ID.
    /// - to      :   Node ID.
    pub fn add_edge(&mut self, from: usize, to: usize) {
        assert!(self.nodes.contains(&from), "A node with ID {} does not exist.", from);
        assert!(self.nodes.contains(&to), "A node with ID {} does not exist.", to);
        self.outgoing_edges.entry(from).or_insert_with(HashSet::new).insert(to);
        self.outgoing_edges.entry(to).or_insert_with(HashSet::new).insert(from);
        self.incoming_edges.entry(from).or_insert_with(HashSet::new).insert(to);
        self.incoming_edges.entry(to).or_insert_with(HashSet::new).insert(from);
    }

    /// Adds an attribute to an edge.
    ///
    /// Parameters:
    ///     from    :   Node ID.
    ///     to      :   Node ID.
    ///     key     :   Attribute key.
    ///     value   :   Attribute value.
    pub fn add_edge_attribute(&mut self, from: usize, to: usize, key: &str, value: Box<dyn Any>) {
        match self.outgoing_edges.get(&from) {
            Some(destinations) => {
                assert!(destinations.contains(&to), "The edge {}-->{} does not exist.", from, to);
                self.edge_attributes
                    .entry((from, to))
                    .or_insert_with(HashMap::new)
                    .insert(key.into(), clone_boxed_any(&value));
            }
            None => {
                panic!("The edge {}-->{} does not exist.", from, to);
            }
        }
    }

    /// Add a new node.
    ///
    /// Parameters:
    /// - id      :   Node ID.
    ///
    /// Returns:
    /// - The new node ID.
    pub fn add_new_node(&mut self) -> usize {
        let new_node_id: usize = self.node_id_counter;
        self.node_id_counter += 1;
        self.nodes.insert(new_node_id);
        self.outgoing_edges.entry(new_node_id).or_insert_with(HashSet::new);
        self.incoming_edges.entry(new_node_id).or_insert_with(HashSet::new);
        self.node_attributes.entry(new_node_id).or_insert_with(HashMap::new);
        new_node_id
    }

    /// Adds an attribute to a node.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///     key     :   Attribute key.
    ///     value   :   Attribute value.
    pub fn add_node_attribute(&mut self, id: usize, key: &str, value: Box<dyn Any>) {
        assert!(self.nodes.contains(&id), "A node with ID {} does not exist.", id);
        self.node_attributes
            .entry(id)
            .or_insert_with(HashMap::new)
            .insert(key.into(), clone_boxed_any(&value));
    }

    /// Fetches the attribute value for an edge.
    ///
    /// Parameters:
    ///     from        :   Node ID.
    ///     to          :   Node ID.
    ///     key         :   Attribute key.
    ///
    /// Returns:
    ///     Reference to the attribute value.
    pub fn get_edge_attribute(&self, from: usize, to: usize, key: &str) -> Option<&Box<dyn Any>> {
        self.edge_attributes.get(&(from, to)).and_then(|inner_map| inner_map.get(key))
    }

    /// Fetches the total number of edges.
    ///
    /// Returns:
    ///     Total number of edges.
    pub fn get_edges_count(&self) -> usize {
        let mut unique_edges: HashSet<(usize,usize)> = HashSet::new();
        for from_node_id in self.outgoing_edges.keys() {
            for to_node_id in self.outgoing_edges.get(from_node_id).unwrap() {
                unique_edges.insert((*from_node_id, *to_node_id));
            }
        }
        unique_edges.len()
    }

    /// Fetches a node's incoming degree.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///
    /// Returns:
    ///     Incoming degree.
    pub fn get_in_degree(&self, id: usize) -> usize {
        assert!(self.nodes.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.incoming_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.incoming_edges.get(&id).unwrap().len()
    }

    /// Fetches the attribute value for a node.
    ///
    /// Parameters:
    ///     id          :   Node ID.
    ///     key         :   Attribute key.
    ///
    /// Returns:
    ///     Reference to the attribute value.
    pub fn get_node_attribute(&self, id: usize, key: &str) -> Option<&Box<dyn Any>> {
        assert!(self.nodes.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.node_attributes.contains_key(&id), "A node with ID {} does not exist.", id);
        self.node_attributes.get(&id).and_then(|inner_map| inner_map.get(key))
    }

    /// Get all node IDs.
    ///
    pub fn get_node_ids(&self) -> Vec<usize> {
        self.nodes.iter().cloned().collect()
    }

    /// Fetches the total number of nodes.
    ///
    /// Returns:
    ///     Total node count.
    pub fn get_nodes_count(&self) -> usize {
        self.nodes.len()
    }

    /// Fetches a node's outgoing degree.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///
    /// Returns:
    ///     Outgoing degree.
    pub fn get_out_degree(&self, id: usize) -> usize {
        assert!(self.nodes.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.outgoing_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.outgoing_edges.get(&id).unwrap().len()
    }

    /// Removes a node.
    ///
    /// Parameters:
    ///     id    :   Node ID.
    pub fn remove_node(&mut self, id: usize) {
        self.nodes.remove(&id);
        self.outgoing_edges.remove(&id);
        self.incoming_edges.remove(&id);
        for (_, set) in &mut self.outgoing_edges {
            set.retain(|x| *x != id);
        }
        for (_, set) in &mut self.incoming_edges {
            set.retain(|x| *x != id);
        }
        self.node_attributes.remove(&id);
        self.edge_attributes.retain(|(first, second), _| *first != id && *second != id);
    }
}
