use exacto_core::prelude::Strand;
use std::collections::HashSet;

use super::*;


#[test]
fn variant_call_returns_consensus_variant_record() {
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGACT".into(),
        VariantType::Insertion
    );
    
    let variant_record_1: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_1.clone()
    );
    
    let variant_record_2: VariantRecord = VariantRecord::new(
        2,
        0,
        1,
        go_1.clone()
    );

    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    variant_records.insert(variant_record_1);
    variant_records.insert(variant_record_2);
    
    let mut variant_call: VariantCall = VariantCall::from_variant_records(
        1,
        variant_records,
        0,
        4,
        6,
        2
    );

    assert_eq!(variant_call.get_consensus_read_ids().len(), 2);
}


/// A read with two records in a call is one read of the call.
///
///   Read   Records
///   1      deletions at 1,000 and at 1,002
///   2      deletion at 1,000
///   3      deletion at 1,000
#[test]
fn get_num_reads_returns_number_of_reads() {
    let deletion = |position: u32| -> GraphOperation {
        GraphOperation::new(
            0,
            position - 1,
            Strand::Forward,
            GraphOperationType::Downstream,
            0,
            position + 1,
            Strand::Forward,
            GraphOperationType::Upstream,
            "".into(),
            VariantType::Deletion
        )
    };
    let variant_records: HashSet<VariantRecord> = HashSet::from([
        VariantRecord::new(1, 100, 101, deletion(1_000)),
        VariantRecord::new(1, 101, 102, deletion(1_002)),
        VariantRecord::new(2, 80, 81, deletion(1_000)),
        VariantRecord::new(3, 90, 91, deletion(1_000))
    ]);

    let variant_call: VariantCall = VariantCall::from_variant_records(
        1,
        variant_records,
        0,
        4,
        6,
        2
    );

    // Read 1 holds two records of the call and is one read (F-12 of the code review).
    assert_eq!(variant_call.get_variant_records().len(), 4);
    assert_eq!(variant_call.get_read_ids(), vec![1, 2, 3]);
    assert_eq!(variant_call.get_num_reads(), 3);
}


/// A constructed consensus is stored in the orientation of its strand.
///
/// No two records spell the insertion the same way, so the consensus is constructed.
/// All three reads are on the reverse strand, so each record holds the reverse complement
/// of the spelling below.
///
///   Read   Spelling in forward orientation
///   1      AAAACCCCGGGGTTTTACGTAC
///   2      AAAACCCCGGGGTTTTACGTAG
///   3      AAAACCCCGGGCTTTTACGTAC
#[test]
fn variant_call_returns_constructed_consensus_in_forward_orientation_for_reverse_strand_records() {
    let variant_records: HashSet<VariantRecord> = HashSet::from([
        VariantRecord::new(1, 500, 521, GraphOperation::new(
            0,
            7_674_700,
            Strand::Reverse,
            GraphOperationType::Downstream,
            0,
            7_674_701,
            Strand::Reverse,
            GraphOperationType::Upstream,
            "GTACGTAAAACCCCGGGGTTTT".into(),
            VariantType::Insertion
        )),
        VariantRecord::new(2, 500, 521, GraphOperation::new(
            0,
            7_674_700,
            Strand::Reverse,
            GraphOperationType::Downstream,
            0,
            7_674_701,
            Strand::Reverse,
            GraphOperationType::Upstream,
            "CTACGTAAAACCCCGGGGTTTT".into(),
            VariantType::Insertion
        )),
        VariantRecord::new(3, 500, 521, GraphOperation::new(
            0,
            7_674_700,
            Strand::Reverse,
            GraphOperationType::Downstream,
            0,
            7_674_701,
            Strand::Reverse,
            GraphOperationType::Upstream,
            "GTACGTAAAAGCCCGGGGTTTT".into(),
            VariantType::Insertion
        ))
    ]);

    let variant_call: VariantCall = VariantCall::from_variant_records(
        1,
        variant_records,
        0,
        4,
        6,
        2
    );

    assert_eq!(*variant_call.get_consensus_method(), GraphOperationConsensusMethod::Constructed);
    assert_eq!(*variant_call.get_consensus_graph_operation().get_strand_1(), Strand::Reverse);
    assert_eq!(variant_call.get_consensus_graph_operation().get_sequence(), "GTACGTAAAACCCCGGGGTTTT");
    assert_eq!(variant_call.get_consensus_graph_operation().get_standardized_sequence(), "AAAACCCCGGGGTTTTACGTAC");
}


/// A constructed consensus does not depend on the order of a `HashSet`.
///
/// The same three spellings as above, read 1 on the forward strand and reads 2 and 3 on the
/// reverse strand. The call is built 40 times; a `HashSet` of three records has six orders.
#[test]
fn variant_call_returns_same_constructed_consensus_for_mixed_strand_records() {
    for _ in 0..40 {
        let variant_records: HashSet<VariantRecord> = HashSet::from([
            VariantRecord::new(1, 500, 521, GraphOperation::new(
                0,
                7_674_700,
                Strand::Forward,
                GraphOperationType::Downstream,
                0,
                7_674_701,
                Strand::Forward,
                GraphOperationType::Upstream,
                "AAAACCCCGGGGTTTTACGTAC".into(),
                VariantType::Insertion
            )),
            VariantRecord::new(2, 500, 521, GraphOperation::new(
                0,
                7_674_700,
                Strand::Reverse,
                GraphOperationType::Downstream,
                0,
                7_674_701,
                Strand::Reverse,
                GraphOperationType::Upstream,
                "CTACGTAAAACCCCGGGGTTTT".into(),
                VariantType::Insertion
            )),
            VariantRecord::new(3, 500, 521, GraphOperation::new(
                0,
                7_674_700,
                Strand::Reverse,
                GraphOperationType::Downstream,
                0,
                7_674_701,
                Strand::Reverse,
                GraphOperationType::Upstream,
                "GTACGTAAAAGCCCGGGGTTTT".into(),
                VariantType::Insertion
            ))
        ]);

        let variant_call: VariantCall = VariantCall::from_variant_records(
            1,
            variant_records,
            0,
            4,
            6,
            2
        );

        // The strand is that of the record with the lowest read ID.
        assert_eq!(*variant_call.get_consensus_method(), GraphOperationConsensusMethod::Constructed);
        assert_eq!(*variant_call.get_consensus_graph_operation().get_strand_1(), Strand::Forward);
        assert_eq!(variant_call.get_consensus_graph_operation().get_sequence(), "AAAACCCCGGGGTTTTACGTAC");
        assert_eq!(variant_call.get_consensus_graph_operation().get_standardized_sequence(), "AAAACCCCGGGGTTTTACGTAC");
    }
}
