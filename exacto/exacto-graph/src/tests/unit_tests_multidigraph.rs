use std::collections::HashSet;

use crate::structs::multidigraph::MultiDiGraph;


#[test]
fn test_multidigraph_add_new_node_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let result: usize = multidigraph.get_nodes_count();
    assert_eq!(result, 3);
}

#[test]
fn test_multidigraph_add_edge_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    let result: usize = multidigraph.get_edges_count();
    assert_eq!(result, 1);
}

#[test]
fn test_multidigraph_add_edge_attribute_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge_attribute(node_1, node_2, "key", Box::new(1usize));
    let result: usize = *multidigraph.get_edge_attribute(node_1, node_2, "key").unwrap().downcast_ref::<usize>().unwrap();
    assert_eq!(result, 1);
}

#[test]
fn test_multidigraph_add_node_attribute_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_node_attribute(node_1, "key", Box::new(1usize));
    let result: usize = *multidigraph.get_node_attribute(node_1, "key").unwrap().downcast_ref::<usize>().unwrap();
    assert_eq!(result, 1);
}

#[test]
fn test_multidigraph_find_paths_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_1, node_3);
    let paths: Vec<Vec<usize>> = multidigraph.find_paths(node_1, node_4);
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_multidigraph_find_paths_2() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_1);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_2);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_4, node_3);
    multidigraph.add_edge(node_1, node_3);
    multidigraph.add_edge(node_3, node_1);
    let paths: Vec<Vec<usize>> = multidigraph.find_paths(node_1, node_4);
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_multidigraph_find_subgraphs_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    let node_5: usize = multidigraph.add_new_node();
    let node_6: usize = multidigraph.add_new_node();
    let node_7: usize = multidigraph.add_new_node();
    let node_8: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_5, node_6);
    multidigraph.add_edge(node_6, node_7);
    multidigraph.add_edge(node_7, node_8);
    let subgraphs: Vec<(usize, HashSet<usize>)> = multidigraph.find_subgraphs();
    assert_eq!(subgraphs.len(), 2);
}

#[test]
fn test_multidigraph_get_incoming_node_ids_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_1, node_3);
    let incoming_node_ids: HashSet<usize> = multidigraph.get_incoming_node_ids(node_3);
    assert_eq!(incoming_node_ids.len(), 2);
}

#[test]
fn test_multidigraph_get_node_levels_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_1, node_3);
    let node_levels: Vec<(usize,usize,usize)> = multidigraph.get_node_levels();
    assert_eq!(node_levels.len(), 4);
}

#[test]
fn test_multidigraph_get_source_node_ids_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_1, node_3);
    let source_node_ids: Vec<usize> = multidigraph.get_source_node_ids();
    assert_eq!(source_node_ids.len(), 1);
    assert_eq!(source_node_ids[0], node_1);
}

#[test]
fn test_multidigraph_get_sink_node_ids_1() {
    let mut multidigraph: MultiDiGraph = MultiDiGraph::new();
    let node_1: usize = multidigraph.add_new_node();
    let node_2: usize = multidigraph.add_new_node();
    let node_3: usize = multidigraph.add_new_node();
    let node_4: usize = multidigraph.add_new_node();
    multidigraph.add_edge(node_1, node_2);
    multidigraph.add_edge(node_2, node_3);
    multidigraph.add_edge(node_3, node_4);
    multidigraph.add_edge(node_1, node_3);
    let sink_node_ids: Vec<usize> = multidigraph.get_sink_node_ids();
    assert_eq!(sink_node_ids.len(), 1);
    assert_eq!(sink_node_ids[0], node_4);
}
