use exacto_core::prelude::LIST_SEPARATOR;
use std::collections::HashSet;

use crate::prelude::*;

use crate::io::builders::build_consensus_sequence_records;
use crate::io::dataframes::consensus_sequence_records_to_dataframe;


fn consensus_set() -> ConsensusSequenceSet {
    let mut set: ConsensusSequenceSet = ConsensusSequenceSet::new();
    // Insertion order deliberately not sorted, so a builder that emitted the set's own order
    // would be caught.
    set.add(ConsensusSequence::new(
        7,
        "ACGT".into(),
        HashSet::from(["read-c/1/ccs".into(), "read-a/1/ccs".into(), "read-b/1/ccs".into()])
    ));
    set.add(ConsensusSequence::new(0, "TTTT".into(), HashSet::from(["solo/1/ccs".into()])));
    set
}


/// The reads a consensus was called from reach the output table, not just their count.
#[test]
fn test_consensus_records_carry_read_names() {
    let records: Vec<ConsensusSequenceRecord> =
        build_consensus_sequence_records(&consensus_set()).collect();
    assert_eq!(records.len(), 2);

    assert_eq!(records[0].cluster_id, 7);
    assert_eq!(records[0].num_reads, 3);
    assert_eq!(
        records[0].read_names.as_ref(),
        ["read-a/1/ccs", "read-b/1/ccs", "read-c/1/ccs"].join(LIST_SEPARATOR),
        "read names must be emitted sorted, not in HashSet iteration order"
    );

    assert_eq!(records[1].read_names.as_ref(), "solo/1/ccs");

    // num_reads and read_names describe the same set; a join that silently dropped a name
    // would otherwise pass the assertions above on a longer list.
    for record in records.iter() {
        assert_eq!(
            record.read_names.split(LIST_SEPARATOR).count(),
            record.num_reads
        );
    }
}


/// Separately-built sets holding the same reads produce byte-identical rows.
///
/// Each `consensus_set()` call builds fresh `HashSet`s, and every `HashSet` seeds its own
/// hasher, so the same reads iterate in a different order per instance. Re-iterating one
/// `HashSet` would not: its order is fixed once built, which would make this test pass against
/// an implementation that never sorted at all.
#[test]
fn test_consensus_records_are_deterministic() {
    let first: Vec<ConsensusSequenceRecord> =
        build_consensus_sequence_records(&consensus_set()).collect();
    for _ in 0..32 {
        let again: Vec<ConsensusSequenceRecord> =
            build_consensus_sequence_records(&consensus_set()).collect();
        assert_eq!(first, again);
    }
}


/// The dataframe shape follows the record shape — `determine-rna-consensus` can return either.
#[test]
fn test_consensus_dataframe_exposes_read_names() {
    let dataframe = consensus_sequence_records_to_dataframe(
        build_consensus_sequence_records(&consensus_set()).collect::<Vec<_>>()
    );
    assert_eq!(
        dataframe.get_column_names_str(),
        vec!["cluster_id", "consensus_sequence", "num_reads", "read_names"]
    );
    assert_eq!(dataframe.height(), 2);
}
