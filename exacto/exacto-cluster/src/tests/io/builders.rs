use std::collections::{HashMap, HashSet};
use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;

use super::*;


/// A call writes one row, its consensus operation, however many spellings its reads hold.
///
/// Eight reads of one cluster carry a 12-base insertion after chr17:7,674,224, a 30-base
/// deletion of chr17:7,674,301-7,674,330 and a breakpoint joining chr17:7,675,000 to
/// chr17:7,680,000:
///
///   Reads   Insertion            Deletion                Breakpoint
///   0-5     CCCATCCGCCTG         7,674,300-7,674,331     7,675,000-7,680,000
///   6       CCATCCGCCTG          7,674,300-7,674,331     7,675,000-7,680,000
///   7       CCCATCCGCCCTG        7,674,301-7,674,331     7,675,002-7,680,000
///
/// Reads 6 and 7 spell the events with a sequencing error.
#[test]
fn build_rna_read_cluster_variant_records_returns_consensus_operations() {
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for read_id in 0..8usize {
        read_names_map.insert(format!("read_{}", read_id).into(), read_id);
    }
    let insertion = |sequence: &str| -> GraphOperation {
        GraphOperation::new(
            0, 7_674_224, Strand::Reverse, GraphOperationType::Downstream,
            0, 7_674_225, Strand::Reverse, GraphOperationType::Upstream,
            sequence.into(), VariantType::Insertion
        )
    };
    let deletion = |position_1: u32| -> GraphOperation {
        GraphOperation::new(
            0, position_1, Strand::Reverse, GraphOperationType::Downstream,
            0, 7_674_331, Strand::Reverse, GraphOperationType::Upstream,
            "".into(), VariantType::Deletion
        )
    };
    let breakpoint = |position_1: u32| -> GraphOperation {
        GraphOperation::new(
            0, position_1, Strand::Reverse, GraphOperationType::Downstream,
            0, 7_680_000, Strand::Reverse, GraphOperationType::Upstream,
            "".into(), VariantType::Breakpoint
        )
    };
    let mut insertion_records: HashSet<VariantRecord> = HashSet::new();
    let mut deletion_records: HashSet<VariantRecord> = HashSet::new();
    let mut breakpoint_records: HashSet<VariantRecord> = HashSet::new();
    for read_id in 0..7usize {
        breakpoint_records.insert(VariantRecord::new(read_id, 300, 301, breakpoint(7_675_000)));
    }
    breakpoint_records.insert(VariantRecord::new(7, 300, 301, breakpoint(7_675_002)));
    for read_id in 0..6usize {
        insertion_records.insert(VariantRecord::new(read_id, 100, 113, insertion("CCCATCCGCCTG")));
        deletion_records.insert(VariantRecord::new(read_id, 200, 201, deletion(7_674_300)));
    }
    insertion_records.insert(VariantRecord::new(6, 100, 112, insertion("CCATCCGCCTG")));
    insertion_records.insert(VariantRecord::new(7, 100, 114, insertion("CCCATCCGCCCTG")));
    deletion_records.insert(VariantRecord::new(6, 200, 201, deletion(7_674_300)));
    deletion_records.insert(VariantRecord::new(7, 200, 201, deletion(7_674_301)));

    let mut cluster_set: RNAReadClusterSet = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        1,
        (0..8usize).collect(),
        Vec::new(),
        vec![
            VariantCall::from_variant_records(0, insertion_records, 0, 4, 6, 2),
            VariantCall::from_variant_records(1, deletion_records, 0, 4, 6, 2),
            VariantCall::from_variant_records(2, breakpoint_records, 0, 4, 6, 2)
        ],
        HashMap::new(),
        HashSet::new()
    ));

    let rows: Vec<(usize, u32, u32, Box<str>, Box<str>)> = build_rna_read_cluster_variant_records(&cluster_set)
        .map(|record| (record.cluster_id, record.position_1, record.position_2, record.sequence, record.variant_type))
        .collect();

    assert_eq!(
        rows,
        vec![
            (1, 7_674_224, 7_674_225, "CCCATCCGCCTG".into(), "INS".into()),
            (1, 7_674_300, 7_674_331, "".into(), "DEL".into()),
            (1, 7_675_000, 7_680_000, "".into(), "BND".into())
        ]
    );
}
