use std::collections::HashSet;

use crate::graph::multidigraph::MultiDiGraph;


#[test]
fn test_multidigraph_add_new_node_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let result: usize = multidigraph.get_nodes_count();
    assert_eq!(result, 3);
}

#[test]
fn test_multidigraph_add_edge_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    multidigraph.add_edge(node_1, node_2,1);
    let result: usize = multidigraph.get_edges_count();
    assert_eq!(result, 1);
}

#[test]
fn test_multidigraph_add_edge_attribute_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    multidigraph.add_edge(node_1, node_2, 100);
    let result: u32 = *multidigraph.get_edge_data(node_1, node_2);
    assert_eq!(result, 100);
}

#[test]
fn test_multidigraph_add_node_attribute_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(100);
    let node_2: usize = multidigraph.add_node(200);
    multidigraph.add_edge(node_1, node_2, 1000);
    let result: u32 = *multidigraph.get_node_data(node_1);
    assert_eq!(result, 100);
}

#[test]
fn test_multidigraph_find_paths_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_1, node_3, 4);
    let paths: Vec<Vec<usize>> = multidigraph.find_paths(node_1, node_4);
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_multidigraph_find_paths_2() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_1, 2);
    multidigraph.add_edge(node_2, node_3, 3);
    multidigraph.add_edge(node_3, node_2, 4);
    multidigraph.add_edge(node_3, node_4, 5);
    multidigraph.add_edge(node_4, node_3, 6);
    multidigraph.add_edge(node_1, node_3, 7);
    multidigraph.add_edge(node_3, node_1, 8);
    let paths: Vec<Vec<usize>> = multidigraph.find_paths(node_1, node_4);
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_multidigraph_find_subgraphs_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    let node_5: usize = multidigraph.add_node(5);
    let node_6: usize = multidigraph.add_node(6);
    let node_7: usize = multidigraph.add_node(7);
    let node_8: usize = multidigraph.add_node(8);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_5, node_6, 4);
    multidigraph.add_edge(node_6, node_7, 5);
    multidigraph.add_edge(node_7, node_8, 6);
    let subgraphs: Vec<(usize, HashSet<usize>)> = multidigraph.find_subgraphs();
    assert_eq!(subgraphs.len(), 2);
}

#[test]
fn test_multidigraph_get_incoming_node_ids_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_1, node_3, 4);
    let incoming_node_ids: &HashSet<usize> = multidigraph.get_incoming_node_ids(node_3);
    assert_eq!(incoming_node_ids.len(), 2);
}

#[test]
fn test_multidigraph_get_node_levels_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_1, node_3, 4);
    let node_levels: Vec<(usize,usize,usize)> = multidigraph.get_node_levels();
    assert_eq!(node_levels.len(), 4);
}

#[test]
fn test_multidigraph_get_source_node_ids_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_1, node_3, 4);
    let source_node_ids: Vec<usize> = multidigraph.get_source_node_ids();
    assert_eq!(source_node_ids.len(), 1);
    assert_eq!(source_node_ids[0], node_1);
}

#[test]
fn test_multidigraph_get_sink_node_ids_1() {
    let mut multidigraph: MultiDiGraph<u32, u32> = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_node(1);
    let node_2: usize = multidigraph.add_node(2);
    let node_3: usize = multidigraph.add_node(3);
    let node_4: usize = multidigraph.add_node(4);
    multidigraph.add_edge(node_1, node_2, 1);
    multidigraph.add_edge(node_2, node_3, 2);
    multidigraph.add_edge(node_3, node_4, 3);
    multidigraph.add_edge(node_1, node_3, 4);
    let sink_node_ids: Vec<usize> = multidigraph.get_sink_node_ids();
    assert_eq!(sink_node_ids.len(), 1);
    assert_eq!(sink_node_ids[0], node_4);
}
