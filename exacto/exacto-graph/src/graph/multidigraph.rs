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
use std::collections::{HashMap,HashSet,VecDeque};

use crate::prelude::*;


/// MultiDiGraph
#[derive(Clone,Debug)]
pub struct MultiDiGraph<N, E> {
    /// Node IDs in the graph.
    pub node_ids: HashSet<usize>,

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

    /// Node data.
    ///
    /// HashMap
    ///     Key     :   Node ID.
    ///     Value   :   Node data.
    pub node_data: HashMap<usize, N>,

    /// Edge data.
    /// 
    /// HashMap
    ///     Key     :   Tuple (source node ID, destination node ID).
    ///     Value   :   Edge data.
    pub edge_data: HashMap<(usize,usize), E>,
    
    pub next_node_id: usize
}

impl<N, E> MultiDiGraph<N, E> {
    /// Constructor
    pub fn new() -> Self {
        MultiDiGraph {
            node_ids: HashSet::new(),
            outgoing_edges: HashMap::new(),
            incoming_edges: HashMap::new(),
            node_data: HashMap::new(),
            edge_data: HashMap::new(),
            next_node_id: 1
        }
    }

    /// Add a new node.
    ///
    /// Parameters:
    /// - id      :   Node ID.
    ///
    /// Returns:
    /// - The new node ID.
    pub fn add_node(&mut self, data: N) -> usize {
        let new_node_id: usize = self.next_node_id;
        self.next_node_id += 1;
        self.node_ids.insert(new_node_id);
        self.outgoing_edges.entry(new_node_id).or_insert_with(HashSet::new);
        self.incoming_edges.entry(new_node_id).or_insert_with(HashSet::new);
        self.node_data.insert(new_node_id, data);
        new_node_id
    }
    
    /// Adds an edge.
    ///
    /// Parameters:
    ///     from        :   Node ID.
    ///     to          :   Node ID.
    pub fn add_edge(&mut self, from: usize, to: usize, data: E) {
        assert_eq!(self.node_ids.contains(&from), true, "A node with ID {} does not exist.", from);
        assert_eq!(self.node_ids.contains(&to), true, "A node with ID {} does not exist.", to);
        self.outgoing_edges.entry(from).or_insert_with(HashSet::new).insert(to);
        self.incoming_edges.entry(to).or_insert_with(HashSet::new).insert(from);
        self.edge_data.insert((from,to), data);
    }

    /// Get layers by breadth-first searach.
    ///
    /// Parameters:
    /// - source    :   Source node ID.
    /// - direction :   "outgoing" or "incoming"
    pub fn bfs_layers(&self, source: usize, direction: MultiDiGraphDirections) -> HashMap<usize, usize> {
        let mut levels: HashMap<usize,usize> = HashMap::new();
        let mut visited: HashSet<usize> = HashSet::new();
        let mut queue: VecDeque<usize> = VecDeque::new();

        // Initialize BFS with the source node
        visited.insert(source);
        queue.push_back(source);
        levels.insert(source, 1); // Source node is at level 1

        // Perform BFS
        while let Some(node_id) = queue.pop_front() {
            let current_level = *levels.get(&node_id).unwrap();

            // Traverse all neighbors of the current node
            let neighbor_node_ids: &HashSet<usize>;
            if direction == MultiDiGraphDirections::Outgoing {
                neighbor_node_ids = self.get_outgoing_node_ids(node_id);
            } else {
                neighbor_node_ids = self.get_incoming_node_ids(node_id);
            }
            for neighbor in neighbor_node_ids.iter() {
                // If the neighbor hasn't been visited yet
                if !visited.contains(neighbor) {
                    visited.insert(neighbor.clone());
                    queue.push_back(neighbor.clone());
                    levels.insert(neighbor.clone(), current_level + 1);
                }
            }
        }

        levels
    }

    /// Displays the MultiDiGraph.
    pub fn display(&self) {
        println!("-------------------------------------------");
        println!("MultiDiGraph ({} nodes and {} edges)", self.node_ids.len(), self.outgoing_edges.keys().len());
        // println!("Is acyclic: {}", self.has_eulerian_cycle());
        println!("-------------------------------------------");
        println!("Node IDs");
        println!("-------------------------------------------");
        for node_id in self.node_ids.iter() {
            println!("{}", node_id);
        }
        println!("-------------------------------------------");
        println!("From                 | To                  ");
        println!("-------------------------------------------");
        for (from_node_id,destinations) in self.outgoing_edges.iter() {
            for to_node_id in destinations.iter() {
                println!("{:20} | {:20}", from_node_id, to_node_id);
            }
        }
        println!("-------------------------------------------");
    }
    
    pub fn edge_exists(&self, from: usize, to: usize) -> bool {
        self.outgoing_edges.contains_key(&from) && self.outgoing_edges[&from].contains(&to)
    }

    /// Finds all paths between two nodes.
    ///
    /// Parameters:
    ///     start     :   Node ID.
    ///     end       :   Node ID.
    ///
    /// Returns:
    ///     A vector of paths where each path is a vector of node IDs.
    pub fn find_paths(&self, start: usize, end: usize) -> Vec<Vec<usize>> {
        let mut all_paths: Vec<Vec<usize>> = Vec::new();
        let mut current_path: Vec<usize> = Vec::new();
        let mut visited: HashSet<usize> = HashSet::new();
        current_path.push(start);
        visited.insert(start);
        self.find_paths_dfs(
            start,
            end,
            &mut visited,
            &mut current_path,
            &mut all_paths
        );
        all_paths
    }

    /// Finds all subgraphs.
    ///
    /// Returns:
    ///     A vector of tuples where each tuple is (subgraph ID, a vector of node IDs).
    pub fn find_subgraphs(&self) -> Vec<(usize, HashSet<usize>)> {
        // Step 1. Identify clusters of nodes
        let mut uf: UnionFind = UnionFind::new();
        for node_id in self.node_ids.iter() {
            uf.union(*node_id as u32, *node_id as u32);
            for &neighbor in self.get_outgoing_node_ids(*node_id) {
                uf.union(*node_id as u32, neighbor as u32);
            }
            for &neighbor in self.get_incoming_node_ids(*node_id) {
                uf.union(*node_id as u32, neighbor as u32);
            }
        }

        // Step 2. Identify subgraphs
        let mut subgraph_idx: usize = 1;
        let mut subgraphs: Vec<(usize, HashSet<usize>)> = Vec::new();
        for cluster in uf.get_clusters().iter() {
            subgraphs.push((subgraph_idx, cluster.iter().map(|&x| x as usize).collect()));
            subgraph_idx += 1;
        }

        subgraphs
    }

    /// Fetches the data for an edge.
    ///
    /// Parameters:
    ///     from        :   Node ID.
    ///     to          :   Node ID.
    ///
    /// Returns:
    ///     Reference to the attribute value.
    pub fn get_edge_data(&self, from: usize, to: usize) -> &E {
        self.edge_data.get(&(from, to)).unwrap()
    }

    /// Fetches the data for an edge.
    ///
    /// Parameters:
    ///     from        :   Node ID.
    ///     to          :   Node ID.
    ///
    /// Returns:
    ///     Reference to the attribute value.
    pub fn get_edge_data_mut(&mut self, from: usize, to: usize) -> &mut E {
        self.edge_data.get_mut(&(from, to)).unwrap()
    }

    /// Fetches the total number of edges.
    ///
    /// Returns:
    ///     Total number of edges.
    pub fn get_edges_count(&self) -> usize {
        let mut unique_edges: HashSet<(usize, usize)> = HashSet::new();
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
        assert!(self.node_ids.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.incoming_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.incoming_edges.get(&id).unwrap().len()
    }

    /// Fetches a node's incoming node IDs.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///
    /// Returns:
    ///     HashSet of incoming node IDs.
    pub fn get_incoming_node_ids(&self, id: usize) -> &HashSet<usize> {
        assert!(self.node_ids.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.incoming_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.incoming_edges.get(&id).as_ref().unwrap()
    }

    /// Fetches the total number of nodes in this MultiDiGraph.
    ///
    /// Returns:
    ///     Total number of nodes.
    pub fn get_nodes_count(&self) -> usize {
        self.node_ids.len()
    }

    /// Fetches the data for a node.
    ///
    /// Parameters:
    ///     node_id     :   Node ID.
    ///
    /// Returns:
    ///     Reference to the node data.
    pub fn get_node_data(&self, node_id: usize) -> &N {
        self.node_data.get(&node_id).unwrap()
    }

    /// Fetches the mutable data for a node.
    ///
    /// Parameters:
    ///     node_id     :   Node ID.
    ///
    /// Returns:
    ///     Reference to the mutable node data.
    pub fn get_node_data_mut(&mut self, node_id: usize) -> &mut N {
        self.node_data.get_mut(&node_id).unwrap()
    }

    /// Get all node IDs.
    ///
    pub fn get_node_ids(&self) -> Vec<usize> {
        self.node_ids.iter().cloned().collect()
    }
    
    /// Fetches node levels.
    ///
    /// Returns:
    ///     Vector of tuples where each tuple is (subgraph ID, node level, node ID).
    pub fn get_node_levels(&self) -> Vec<(usize,usize,usize)> {
        // Step 1. Find all subgraphs
        let subgraphs: Vec<(usize, HashSet<usize>)> = self.find_subgraphs();

        // Step 2. Get level data
        let mut node_levels: Vec<(usize,usize,usize)> = Vec::new();
        let source_node_ids: HashSet<usize> = self.get_source_node_ids().into_iter().collect();
        for (subgraph_id, subgraph_node_ids) in subgraphs.iter() {
            // Prepare a mutable HashMap of the nodes' incoming degrees
            let mut in_degree: HashMap<usize,usize> = HashMap::new();
            for subgraph_node_id in subgraph_node_ids.iter() {
                in_degree.insert(*subgraph_node_id, self.get_incoming_node_ids(*subgraph_node_id).len());
            }

            // Append all source nodes into the queue
            let mut queue = VecDeque::new();
            for subgraph_node_id in subgraph_node_ids.iter() {
                if source_node_ids.contains(subgraph_node_id) {
                    queue.push_back((subgraph_node_id.clone(), 1));
                }
            }

            // Compute the node level by Kahn's algorithm (topological sorting)
            while let Some((node_id,node_level)) = queue.pop_front() {
                node_levels.push((*subgraph_id, node_level, node_id));
                let outgoing_node_ids: &HashSet<usize> = self.get_outgoing_node_ids(node_id);
                for outgoing_node_id in outgoing_node_ids.iter() {
                    if let Some(degree) = in_degree.get_mut(outgoing_node_id) {
                        *degree -= 1;
                        if *degree == 0 {
                            queue.push_back((*outgoing_node_id, node_level + 1));
                        }
                    }
                }
            }
        }

        node_levels
    }

    /// Fetches a node's outgoing degree.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///
    /// Returns:
    ///     Outgoing degree.
    pub fn get_out_degree(&self, id: usize) -> usize {
        assert!(self.node_ids.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.outgoing_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.outgoing_edges.get(&id).unwrap().len()
    }

    /// Gets a node's outgoing node IDs.
    ///
    /// Parameters:
    ///     id      :   Node ID.
    ///
    /// Returns:
    ///     HashSet of outgoing node IDs.
    pub fn get_outgoing_node_ids(&self, id: usize) -> &HashSet<usize> {
        assert!(self.node_ids.contains(&id), "A node with ID {} does not exist.", id);
        assert!(self.outgoing_edges.contains_key(&id), "A node with ID {} does not exist.", id);
        self.outgoing_edges.get(&id).as_ref().unwrap()
    }

    /// Fetches all source node IDs.
    ///
    /// Returns:
    ///     Vector of node IDs.
    pub fn get_source_node_ids(&self) -> Vec<usize> {
        let mut source_node_ids: Vec<usize> = Vec::new();
        for node_id in self.node_ids.iter() {
            if self.get_in_degree(*node_id) == 0 {
                source_node_ids.push(*node_id);
            }
        }
        source_node_ids
    }

    /// Fetches all sink node IDs.
    ///
    /// Returns:
    ///     Vector of node IDs.
    pub fn get_sink_node_ids(&self) -> Vec<usize> {
        let mut sink_node_ids: Vec<usize> = Vec::new();
        for node_id in self.node_ids.iter() {
            if self.get_out_degree(*node_id) == 0 {
                sink_node_ids.push(*node_id);
            }
        }
        sink_node_ids
    }
    
    pub fn node_exists(&self, id: usize) -> bool {
        self.node_ids.contains(&id)
    }

    pub fn remove_edge(&mut self, from: usize, to: usize) {
        self.outgoing_edges.get_mut(&from).unwrap().remove(&to);
        self.incoming_edges.get_mut(&to).unwrap().remove(&from);
        self.edge_data.remove(&(from, to));
    }

    /// Removes a node.
    ///
    /// Parameters:
    ///     id    :   Node ID.
    pub fn remove_node(&mut self, id: usize) {
        self.node_ids.remove(&id);
        self.outgoing_edges.remove(&id);
        self.incoming_edges.remove(&id);
        for (_, set) in &mut self.outgoing_edges {
            set.retain(|x| *x != id);
        }
        for (_, set) in &mut self.incoming_edges {
            set.retain(|x| *x != id);
        }
        self.node_data.remove(&id);
        self.edge_data.retain(|(first, second), _| *first != id && *second != id);
    }
}

/// Helper methods
impl<N, E> MultiDiGraph<N, E> {
    /// Finds all paths using depth-first search.
    ///
    /// Parameters:
    ///     current         :   Node ID.
    ///     end             :   Node ID.
    ///     visited         :   HashSet of visited node ID.
    ///     current_path    :   Vector of node IDs.
    ///     all_paths       :   Vector of vector of node IDs.
    fn find_paths_dfs(
        &self,
        current: usize,
        end: usize,
        visited: &mut HashSet<usize>,
        current_path: &mut Vec<usize>,
        all_paths: &mut Vec<Vec<usize>>
    ) {
        if current == end {
            all_paths.push(current_path.clone());
        } else {
            if let Some(neighbors) = self.outgoing_edges.get(&current) {
                for neighbor in neighbors {
                    if !visited.contains(neighbor) {
                        visited.insert(neighbor.clone());
                        current_path.push(neighbor.clone());
                        self.find_paths_dfs(*neighbor, end, visited, current_path, all_paths);
                        current_path.pop();
                        visited.remove(neighbor);
                    }
                }
            }
        }
    }
}
