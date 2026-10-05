use exacto_core::prelude::*;

use super::*;


fn local_indel(position_1: u32, position_2: u32, sequence: &str) -> GraphOperationView {
    GraphOperationView::new(
        "chr1",
        position_1,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr1",
        position_2,
        GraphOperationType::Upstream,
        Strand::Forward,
        sequence,
        None
    )
}

#[test]
fn graph_operation_view_net_length_delta_single_base_deletion_is_negative_one() {
    let gov = local_indel(100, 102, "");
    assert_eq!(gov.net_length_delta(), Some(-1));
    assert!(gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_three_base_deletion_is_inframe() {
    let gov = local_indel(100, 104, "");
    assert_eq!(gov.net_length_delta(), Some(-3));
    assert!(!gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_single_base_insertion_is_positive_one() {
    let gov = local_indel(100, 101, "A");
    assert_eq!(gov.net_length_delta(), Some(1));
    assert!(gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_three_base_insertion_is_inframe() {
    let gov = local_indel(100, 101, "ACG");
    assert_eq!(gov.net_length_delta(), Some(3));
    assert!(!gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_snv_is_zero() {
    let gov = local_indel(100, 102, "A");
    assert_eq!(gov.net_length_delta(), Some(0));
    assert!(!gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_mnv_is_zero() {
    // MNV: gap == sequence length -> balanced substitution, net 0.
    let gov = local_indel(100, 104, "ACG");
    assert_eq!(gov.net_length_delta(), Some(0));
    assert!(!gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_translocation_is_none() {
    let gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr2",
        200,
        GraphOperationType::Upstream,
        Strand::Forward,
        "",
        None
    );
    assert_eq!(gov.net_length_delta(), None);
    assert!(!gov.is_frameshift());
}

#[test]
fn graph_operation_view_net_length_delta_breakpoint_is_none() {
    let gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Upstream,
        Strand::Forward,
        "chr1",
        500,
        GraphOperationType::Downstream,
        Strand::Forward,
        "",
        Some(1)
    );
    assert_eq!(gov.net_length_delta(), None);
    assert!(!gov.is_frameshift());
}
