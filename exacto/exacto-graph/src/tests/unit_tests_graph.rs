use std::any::Any;
use crate::structs::graph::Graph;


#[test]
fn test_graph_add_new_node_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    let result: usize = graph.get_nodes_count();
    assert_eq!(result, 3);
}

#[test]
fn test_graph_add_edge_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_2, node_3);
    let result: usize = graph.get_nodes_count();
    assert_eq!(result, 3);
}

#[test]
fn test_graph_add_edge_2() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    let node_4: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_1, node_3);
    graph.add_edge(node_1, node_4);
    let num_edges: usize = graph.get_out_degree(node_1);
    let num_nodes: usize = graph.get_nodes_count();
    assert_eq!(num_nodes, 4);
    assert_eq!(num_edges, 3);
}

#[test]
fn test_graph_add_node_attribute_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_node_attribute(node_1, "sample_id", Box::new("sample001".to_string().into_boxed_str()));
    let value: Option<&Box<dyn Any>> = graph.get_node_attribute(node_1, "sample_id");
    if let Some(str_value) = value.unwrap().downcast_ref::<Box<str>>() {
        assert_eq!(&**str_value, "sample001");
    } else {
        panic!("The node attribute should be a Box<str> value.");
    }
}

#[test]
fn test_graph_get_node_attribute_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_2, node_3);
    let value: Option<&Box<dyn Any>> = graph.get_node_attribute(node_1, "key");
    assert!(value.is_none(), "The Option should be None, but it was Some.");
}

#[test]
fn test_graph_add_edge_attribute_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge_attribute(node_1,node_2,"key", Box::new(100));
    let value: Option<&Box<dyn Any>> = graph.get_edge_attribute(node_1, node_2, "key");
    if let Some(int_value) = value.unwrap().downcast_ref::<i32>() {
        assert_eq!(*int_value, 100);
    } else {
        panic!("The node attribute should be an integer value.");
    }
}

#[test]
fn test_graph_get_edge_attribute_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    let value: Option<&Box<dyn Any>> = graph.get_edge_attribute(node_1, node_2, "key");
    assert!(value.is_none(), "The Option should be None, but it was Some.");
}

#[test]
fn test_graph_get_in_degree_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    let node_4: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_2, node_3);
    graph.add_edge(node_3, node_4);
    graph.add_edge(node_1, node_3);
    let in_degree: usize = graph.get_in_degree(node_3);
    assert_eq!(in_degree, 3);
}

#[test]
fn test_graph_get_out_degree_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    let node_4: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_2, node_3);
    graph.add_edge(node_3, node_4);
    graph.add_edge(node_1, node_3);
    let in_degree: usize = graph.get_out_degree(node_3);
    assert_eq!(in_degree, 3);
}

#[test]
fn test_graph_remove_node_1() {
    let mut graph: Graph = Graph::new();
    let node_1: usize = graph.add_new_node();
    let node_2: usize = graph.add_new_node();
    let node_3: usize = graph.add_new_node();
    let node_4: usize = graph.add_new_node();
    graph.add_edge(node_1, node_2);
    graph.add_edge(node_2, node_3);
    graph.add_edge(node_3, node_4);
    graph.add_edge(node_1, node_3);
    graph.remove_node(node_2);
    assert_eq!(graph.get_nodes_count(), 3);
    assert_eq!(graph.get_out_degree(node_3), 2);
}
