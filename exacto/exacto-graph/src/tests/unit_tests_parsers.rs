use crate::common::constants::{VarGraphOrientations, VarGraphStrands, VarGraphNodeTypes};
use crate::common::parsers::*;
use crate::structs::vargraph_variant_node::VarGraphVariantNode;
use crate::traits::vargraph_node::VarGraphNode;


#[test]
fn test_parse_graph_operation_1() {
    let variant_node: VarGraphVariantNode = parse_graph_operation("chr1:2:+:D:chr1:3:+:U:AATTG:5");
    assert_eq!(&*variant_node.get_chromosome_1(), "chr1");
    assert_eq!(variant_node.get_position_1(), 2);
    assert_eq!(variant_node.strand_1, VarGraphStrands::Forward);
    assert_eq!(variant_node.orientation_1, VarGraphOrientations::Downstream);
    assert_eq!(&*variant_node.get_chromosome_2(), "chr1");
    assert_eq!(variant_node.get_position_2(), 3);
    assert_eq!(variant_node.strand_2, VarGraphStrands::Forward);
    assert_eq!(variant_node.orientation_2, VarGraphOrientations::Upstream);
    assert_eq!(&*variant_node.get_sequence(), "AATTG");
    assert_eq!(variant_node.get_sequence_length(), 5);
    assert_eq!(variant_node.get_type(), VarGraphNodeTypes::Variant);
}

#[test]
fn test_parse_graph_operation_2() {
    let variant_node: VarGraphVariantNode = parse_graph_operation("chr1:1:+:D:chr1:10:+:U::0");
    assert_eq!(&*variant_node.get_chromosome_1(), "chr1");
    assert_eq!(variant_node.get_position_1(), 1);
    assert_eq!(variant_node.strand_1, VarGraphStrands::Forward);
    assert_eq!(variant_node.orientation_1, VarGraphOrientations::Downstream);
    assert_eq!(&*variant_node.get_chromosome_2(), "chr1");
    assert_eq!(variant_node.get_position_2(), 10);
    assert_eq!(variant_node.strand_2, VarGraphStrands::Forward);
    assert_eq!(variant_node.orientation_2, VarGraphOrientations::Upstream);
    assert_eq!(&*variant_node.get_sequence(), "");
    assert_eq!(variant_node.get_sequence_length(), 0);
    assert_eq!(variant_node.get_type(), VarGraphNodeTypes::Variant);
}

