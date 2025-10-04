use std::fs;
use std::path::Path;
use crate::common::constants::{VarGraphOrientations, VarGraphStrands};
use crate::structs::vargraph::VarGraph;
use crate::structs::vargraph_path::VarGraphPath;
use crate::structs::vargraph_reference_node::VarGraphReferenceNode;


#[test]
fn test_vargraph_from_fasta_file() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    assert_eq!(vargraph.get_nodes_count(), 2);
}

#[test]
fn test_vargraph_add_reference_1() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let reference_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr2", 1, 10, "TTTTCCCCGA"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(vec![&reference_node_1, &reference_node_2]);
    assert_eq!(vargraph.get_nodes_count(), 2);
}

#[test]
fn test_vargraph_add_variant_1() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(vec![&reference_node_1]);
    vargraph.add_variant(
        "chr1",
        4, 
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chr1",
        6,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        "A"
    );
    assert_eq!(vargraph.get_nodes_count(), 6);
}

#[test]
fn test_vargraph_add_variant_2() {
    let reference_node_1: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr1", 1, 10, "ACGTACGTAG"
    );
    let reference_node_2: VarGraphReferenceNode = VarGraphReferenceNode::new(
        "chr2", 1, 10, "CCTGATCGTA"
    );
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(vec![&reference_node_1, &reference_node_2]);
    vargraph.add_variant(
        "chr1",
        5,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chr2",
        6,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        ""
    );
    assert_eq!(vargraph.get_nodes_count(), 7);
}

#[test]
fn test_vargraph_get_linearized_contigs_1() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrA",
        6,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        "T"
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCTTACGTAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_2() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrA",
        11,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        ""
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_3() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrA",
        5,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        "CCC"
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCCCCGTACGTAGCTAGCTAG".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_4() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrB",
        5,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        ""
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCTTCCCAAAGGGTTTCC".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_5() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrB",
        5,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        "CCC"
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCCCCTTCCCAAAGGGTTTCC".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_6() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        15,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Reverse,
        "chrB",
        15,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        ""
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "CTAGCTGTTTCC".into());
}

#[test]
fn test_vargraph_get_linearized_contigs_7() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    vargraph.add_variant(
        "chrA",
        3,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Forward,
        "chrA",
        15,
        VarGraphOrientations::Downstream,
        VarGraphStrands::Reverse,
        ""
    );
    vargraph.add_variant(
        "chrA",
        4,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Reverse,
        "chrA",
        16,
        VarGraphOrientations::Upstream,
        VarGraphStrands::Forward,
        ""
    );
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGTAGCTACGTACGGCTAG".into());
}

