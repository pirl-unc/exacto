use std::io::Write;

use super::*;


/// The records of the RNA consensus TSV file are read, and a file written while the columns
/// `num_iters` and `num_concordant_iters` existed gives the same records.
#[test]
fn load_consensus_sequence_records_returns_records() {
    let temp_dir = tempfile::tempdir().unwrap();

    let tsv_path = temp_dir.path().join("consensus.tsv");
    let mut tsv_file = std::fs::File::create(&tsv_path).unwrap();
    writeln!(tsv_file, "cluster_id\tconsensus_sequence\tnum_reads\tread_names").unwrap();
    writeln!(tsv_file, "0\tTTTT\t1\tsolo/1/ccs").unwrap();
    writeln!(tsv_file, "7\tACGT\t3\tread-a/1/ccs;read-b/1/ccs;read-c/1/ccs").unwrap();
    drop(tsv_file);

    let old_tsv_path = temp_dir.path().join("consensus_old.tsv");
    let mut old_tsv_file = std::fs::File::create(&old_tsv_path).unwrap();
    writeln!(old_tsv_file, "cluster_id\tconsensus_sequence\tnum_reads\tread_names\tnum_iters\tnum_concordant_iters").unwrap();
    writeln!(old_tsv_file, "0\tTTTT\t1\tsolo/1/ccs\t0\t0").unwrap();
    writeln!(old_tsv_file, "7\tACGT\t3\tread-a/1/ccs;read-b/1/ccs;read-c/1/ccs\t5\t3").unwrap();
    drop(old_tsv_file);

    let records: Vec<ConsensusSequenceRecord> = load_consensus_sequence_records(tsv_path.to_str().unwrap());
    let old_records: Vec<ConsensusSequenceRecord> = load_consensus_sequence_records(old_tsv_path.to_str().unwrap());

    assert_eq!(
        records,
        vec![
            ConsensusSequenceRecord {
                cluster_id: 0,
                consensus_sequence: "TTTT".into(),
                num_reads: 1,
                read_names: "solo/1/ccs".into()
            },
            ConsensusSequenceRecord {
                cluster_id: 7,
                consensus_sequence: "ACGT".into(),
                num_reads: 3,
                read_names: "read-a/1/ccs;read-b/1/ccs;read-c/1/ccs".into()
            }
        ]
    );
    assert_eq!(old_records, records);
}
