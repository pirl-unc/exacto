use crate::structs::vargraph::VarGraph;


#[test]
fn test_vargraph_add_reference_1() {
    let mut vargraph: VarGraph = VarGraph::new();
    vargraph.add_reference("chr1",1,10,"ACGTACGTAG");
    vargraph.add_reference("chr2",1,10,"TTTTCCCCGA");
    assert_eq!(vargraph.get_nodes_count(), 2);
}

#[test]
fn test_vargraph_add_variant_1() {
    let mut vargraph: VarGraph = VarGraph::new();
    vargraph.add_reference("chr1",1,10,"ACGTACGTAG");
    vargraph.add_variant("chr1:4:+:D:chr1:6:+:U:A:1");
    assert_eq!(vargraph.get_nodes_count(), 6);
}

#[test]
fn test_vargraph_add_variant_2() {
    let mut vargraph: VarGraph = VarGraph::new();
    vargraph.add_reference("chr1",1,10,"ACGTACGTAG");
    vargraph.add_reference("chr2",1,10,"CCTGATCGTA");
    vargraph.add_variant("chr1:5:+:D:chr2:6:+:U::0");
    assert_eq!(vargraph.get_nodes_count(), 7);
}

