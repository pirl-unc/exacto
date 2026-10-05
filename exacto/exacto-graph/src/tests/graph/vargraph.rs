use std::fs;
use std::path::Path;
use super::*;


#[test]
fn test_vargraph_1() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    assert_eq!(vargraph.get_nodes_count(), 2);
}

#[test]
fn test_vargraph_2() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let reference_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr2", 1, 10, "TTTTCCCCGA"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(
        vec![vec![reference_node_1, reference_node_2]]
    );
    assert_eq!(vargraph.get_nodes_count(), 2);
}

#[test]
fn test_vargraph_3() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(
        vec![vec![reference_node_1]]
    );
    let gov: GraphOperationView = GraphOperationView::new(
        "chr1",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr1",
        6,
        GraphOperationType::Upstream,
        Strand::Forward,
        "A",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    assert_eq!(vargraph.get_nodes_count(), 6);
}

#[test]
fn test_vargraph_4() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let reference_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr2", 1, 10, "CCTGATCGTA"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(
        vec![vec![reference_node_1, reference_node_2]]
    );
    let gov: GraphOperationView = GraphOperationView::new(
        "chr1",
        5,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr2",
        6,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);

    assert_eq!(vargraph.get_nodes_count(), 7);
}

#[test]
fn test_vargraph_5() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrA",
        6,
        GraphOperationType::Upstream,
        Strand::Forward,
        "T",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCTTACGTAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_6() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrA",
        11,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_7() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrA",
        5,
        GraphOperationType::Upstream,
        Strand::Forward,
        "CCC",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCCCCGTACGTAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_8() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrB",
        5,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCTTCCCAAAGGGTTTCC".into());
}

#[test]
fn test_vargraph_9() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrB",
        5,
        GraphOperationType::Upstream,
        Strand::Forward,
        "CCC",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCCCCTTCCCAAAGGGTTTCC".into());
}

#[test]
fn test_vargraph_10() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov: GraphOperationView = GraphOperationView::new(
        "chrA",
        15,
        GraphOperationType::Upstream,
        Strand::Reverse,
        "chrB",
        15,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(1, gov);
    vargraph.add_variant_node(variant_node);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "CTAGCTGTTTCC".into());
}

#[test]
fn test_vargraph_11() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    let gov_1: GraphOperationView = GraphOperationView::new(
        "chrA",
        3,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chrA",
        15,
        GraphOperationType::Downstream,
        Strand::Reverse,
        "",
        None
    );
    let variant_node_1: VarGraphVariantNode = VarGraphVariantNode::new(1, gov_1);
    vargraph.add_variant_node(variant_node_1);
    let gov_2: GraphOperationView = GraphOperationView::new(
        "chrA",
        4,
        GraphOperationType::Upstream,
        Strand::Reverse,
        "chrA",
        16,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    let variant_node_2: VarGraphVariantNode = VarGraphVariantNode::new(2, gov_2);
    vargraph.add_variant_node(variant_node_2);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGTAGCTACGTACGGCTAG".into());
}

#[test]
fn test_vargraph_12() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let reference_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr2", 1, 10, "AGCTAGCTAT"
    );
    let mut vargraph_1: VarGraph = VarGraph::from_reference_nodes(vec![vec![reference_node_1]]);
    let gov_1: GraphOperationView = GraphOperationView::new(
        "chr1",
        4,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr1",
        6,
        GraphOperationType::Upstream,
        Strand::Forward,
        "A",
        None
    );
    let variant_node_1: VarGraphVariantNode = VarGraphVariantNode::new(1, gov_1);
    vargraph_1.add_variant_node(variant_node_1);
    let mut vargraph_2: VarGraph = VarGraph::from_reference_nodes(vec![vec![reference_node_2]]);
    let gov_2: GraphOperationView = GraphOperationView::new(
        "chr2",
        2,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr2",
        4,
        GraphOperationType::Upstream,
        Strand::Forward,
        "T",
        None
    );
    let variant_node_2: VarGraphVariantNode = VarGraphVariantNode::new(2, gov_2);
    vargraph_2.add_variant_node(variant_node_2);
    let mut vargraph: VarGraph = VarGraph::merge(vec![vargraph_1, vargraph_2]);
    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(vargraph.get_nodes_count(), 12);
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_vargraph_split_keeps_disabled_edges() {
    // Review F-07: a deletion chrA:11-29 and an SNV at 20 inside it, added in either order to an
    // unsplit reference node. Splitting the node for the second variant must not re-enable the
    // edges the first disabled, so both orders give the same paths.
    let mut sequences_by_order: Vec<Vec<Box<str>>> = Vec::new();
    for order in [[1usize, 2], [2, 1]] {
        let reference_node: VarGraphReferenceNode = VarGraphReferenceNode::new(
            "chrA", 1, 40, "ACGTACGTACGGGGGGGGGGGGGGGGGGGGTTTTTTTTTT"
        );
        let mut vargraph: VarGraph = VarGraph::from_reference_nodes(vec![vec![reference_node]]);
        for variant_id in order {
            let gov: GraphOperationView = if variant_id == 1 {
                GraphOperationView::new(
                    "chrA", 10, GraphOperationType::Downstream, Strand::Forward,
                    "chrA", 30, GraphOperationType::Upstream, Strand::Forward,
                    "", None
                )
            } else {
                GraphOperationView::new(
                    "chrA", 19, GraphOperationType::Downstream, Strand::Forward,
                    "chrA", 21, GraphOperationType::Upstream, Strand::Forward,
                    "N", None
                )
            };
            vargraph.add_variant_node(VarGraphVariantNode::new(variant_id, gov));
        }
        let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
            &vargraph.get_variant_node_ids().into_iter().collect(),
            &HashSet::new()
        );
        let mut sequences: Vec<Box<str>> = paths.iter().map(|path| path.get_sequence()).collect();
        sequences.sort();
        sequences_by_order.push(sequences);
    }

    assert_eq!(sequences_by_order[0], vec!["ACGTACGTACGTTTTTTTTTT".into(), "GGGGGGGGGNGGGGGGGGG".into()]);
    assert_eq!(sequences_by_order[1], sequences_by_order[0]);
}
