use exacto_core::prelude::Strand;
use std::sync::Arc;

use super::*;



#[test]
fn dna_control_subtraction_filter_returns_false_for_shared_breakpoint_position() {
    // TRA:chr1:1001-chr3:2001
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Downstream,
        2,
        2001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Translocation
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        100,
        100,
        go_1
    );

    // TRA:chr1:995-chr2:1998
    let go_2: GraphOperation = GraphOperation::new(
        0,
        995,
        Strand::Forward,
        GraphOperationType::Downstream,
        1,
        1998,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Translocation
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        500,
        500,
        go_2
    );

    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &[Arc::new(b)],
        100_000,
        0
    );
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        true,
        false
    );

    assert!(control_subtraction_filter.passes(&a) == false);
}

#[test]
fn dna_control_subtraction_filter_returns_false_for_proximal_deletions() {
    // DEL:chr1:1001-1100
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1111,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        100,
        100,
        go_1
    );

    // DEL:chr1:990-1150
    let go_2: GraphOperation = GraphOperation::new(
        0,
        989,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1151,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        300,
        300,
        go_2
    );

    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &[Arc::new(b)],
        100_000,
        0
    );
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        true,
        false
    );

    assert!(control_subtraction_filter.passes(&a) == false);
}

#[test]
fn dna_control_subtraction_filter_returns_false_for_proximal_insertions_1() {
    // INS:chr1:1000
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        31,
        60,
        go_1
    );

    // INS:chr1:1005
    let go_2: GraphOperation = GraphOperation::new(
        0,
        1005,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1006,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CGACTACGATCGACTACGATCGACTACGAT".into(),
        VariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        91,
        120,
        go_2
    );

    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &[Arc::new(b)],
        100_000,
        0
    );
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        true,
        false
    );

    assert!(control_subtraction_filter.passes(&a) == false);
}

#[test]
fn dna_control_subtraction_filter_returns_false_for_proximal_insertions_2() {
    // INS:chr1:1000
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        61,
        90,
        go_1
    );

    // INS:chr1:1002
    let go_2: GraphOperation = GraphOperation::new(
        0,
        1002,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1003,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CGATC".into(),
        VariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        6,
        10,
        go_2
    );

    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &[Arc::new(b)],
        100_000,
        0
    );
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        0.05f64,
        0.95f64,
        1000,
        0.01f64,
        true,
        false
    );

    assert!(control_subtraction_filter.passes(&a) == false);
}

#[test]
fn dna_control_subtraction_filter_returns_true_for_proximal_but_dissimilar_insertions() {
    // INS:chr1:1000
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        31,
        60,
        go_1
    );

    // INS:chr1:1002
    let go_2: GraphOperation = GraphOperation::new(
        0,
        1002,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1003,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        5,
        5,
        go_2
    );

    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &[Arc::new(b)],
        100_000,
        0
    );
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        0.05f64,
        0.95f64,
        1000,
        0.01f64,
        true,
        false
    );

    assert!(control_subtraction_filter.passes(&a) == true);
}
