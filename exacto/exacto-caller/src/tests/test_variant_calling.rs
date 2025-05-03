use bimap::BiMap;
use noodles_bam as bam;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use crate::prelude::*;


#[test]
fn test_variant_call_1() {
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let variant_record_1: VariantRecord = VariantRecord::new(1, go_1.clone());
    let variant_record_2: VariantRecord = VariantRecord::new(2, go_1.clone());
    let mut variant_call: VariantCall = VariantCall::new();
    variant_call.add_variant_record(variant_record_1);
    variant_call.add_variant_record(variant_record_2);

    assert!(variant_call.get_consensus_record().1.len() == 2);
}

#[test]
fn test_cluster_variant_records_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();

    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );

    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_325382_158010/1/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let record: &bam::Record = records_map.get(&read_id).unwrap().first().unwrap();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        read_sequence,
        quality_scores,
        vec![record.clone()]
    );

    let variant_records: Vec<VariantRecord> = alignment.identify_sequence_variant_records(
        30,
        30f32
    );

    assert!(variant_records.len() == 1);

    let variant_calls: Vec<VariantCall> = cluster_variant_records(
        variant_records.iter().map(|record| Arc::new(record.clone())).collect(),
        1,
        0.5f32,
        0.5f32,
        2000,
        1000,
        10000
    );

    assert!(variant_calls.len() == 1);
}

#[test]
fn test_diff_variant_records_1() {
    // INS:chr1:1000
    // DEL:chr1:1100-1200
    // SNV:chr1:1205
    // SNV:chr1:1400
    let mut a: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1001
    // SNV:chr1:1205
    let mut b: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1000
    let go_a1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let vr_a1: VariantRecord = VariantRecord::new(1, go_a1);

    // DEL:chr1:1101-1200
    let go_a2: SequenceOperation = SequenceOperation::new(
        0,
        1100,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let vr_a2: VariantRecord = VariantRecord::new(1, go_a2);

    // SNV:chr1:1205
    let go_a3: SequenceOperation = SequenceOperation::new(
        0,
        1204,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "A".into(),
        SequenceOperationVariantType::SingleNucleotideVariant
    );
    let vr_a3: VariantRecord = VariantRecord::new(1, go_a3);

    // SNV:chr1:1400
    let go_a4: SequenceOperation = SequenceOperation::new(
        0,
        1399,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1401,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "T".into(),
        SequenceOperationVariantType::SingleNucleotideVariant
    );
    let vr_a4: VariantRecord = VariantRecord::new(1, go_a4);

    a.push(vr_a1);
    a.push(vr_a2);
    a.push(vr_a3);
    a.push(vr_a4);

    // INS:chr1:1001
    let go_b1: SequenceOperation = SequenceOperation::new(
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1002,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "CGATCGACTA".into(),
        SequenceOperationVariantType::Insertion
    );
    let vr_b1: VariantRecord = VariantRecord::new(2, go_b1);

    // SNV:chr1:1205
    let go_b2: SequenceOperation = SequenceOperation::new(
        0,
        1204,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "A".into(),
        SequenceOperationVariantType::SingleNucleotideVariant
    );
    let vr_b2: VariantRecord = VariantRecord::new(2, go_b2);

    b.push(vr_b1);
    b.push(vr_b2);

    let a_ref: Vec<Arc<VariantRecord>> = a
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let b_ref: Vec<Arc<VariantRecord>> = b
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let a_only: Vec<Arc<VariantRecord>> = diff_variant_records(
        a_ref,
        b_ref,
        100_000,
        1,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(a_only.len() == 2);

    let mut found_del_1: bool = false;
    let mut found_snv_2: bool = false;

    for variant_record in a_only.iter() {
        if variant_record.get_sequence_operation_boxed_str() == "0:1100:+:D:0:1201:+:U::0:DEL".into() {
            found_del_1 = true;
        }
        if variant_record.get_sequence_operation_boxed_str() == "0:1399:+:D:0:1401:+:U:T:1:SNV".into() {
            found_snv_2 = true;
        }
    }

    assert!(found_del_1);
    assert!(found_snv_2);
}

#[test]
fn test_identify_variants_in_long_read_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();

    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );

    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_202369_785869/3/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        read_sequence,
        quality_scores,
        records_map.get(&read_id).unwrap().clone()
    );

    let variant_records: Vec<VariantRecord> = alignment.identify_sequence_variant_records(
        30,
        30f32
    );

    assert!(variant_records.len() == 1);

    let variant_calls: Vec<VariantCall> = cluster_variant_records(
        variant_records.iter().map(|record| Arc::new(record.clone())).collect(),
        1,
        0.5f32,
        0.5f32,
        2000,
        1000,
        10000
    );

    assert!(variant_calls.len() == 1);
}

#[test]
fn test_is_clusterable_1() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // DEL:chr1:1001-1100
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000
    );

    assert!(result == false);
}

#[test]
fn test_is_clusterable_2() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // INS:chr1:1001
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1002,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "CGATCGACTC".into(),
        SequenceOperationVariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000
    );

    assert!(result == true);
}

#[test]
fn test_is_clusterable_3() {
    // DEL:chr1:1001-1100
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // DEL:chr1:999-1110
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        998,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1111,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000
    );

    assert!(result == true);
}

#[test]
fn test_is_clusterable_4() {
    // TRA:chr1:1001-chr2:2001
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Downstream,
        1,
        2001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Translocation
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // TRA:chr1:995-chr2:1998
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        995,
        Strand::Forward,
        SequenceOperationType::Downstream,
        1,
        1998,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Translocation
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000
    );

    assert!(result == true);
}

#[test]
fn test_is_diffable_1() {
    // TRA:chr1:1001-chr3:2001
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Downstream,
        2,
        2001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Translocation
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // TRA:chr1:995-chr2:1998
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        995,
        Strand::Forward,
        SequenceOperationType::Downstream,
        1,
        1998,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Translocation
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_different(
        &a.sequence_operation,
        &b.sequence_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(result == false);
}

#[test]
fn test_is_diffable_2() {
    // DEL:chr1:1001-1100
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1111,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // DEL:chr1:990-1150
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        989,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1151,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_different(
        &a.sequence_operation,
        &b.sequence_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(result == false);
}

#[test]
fn test_is_diffable_3() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // INS:chr1:1005
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        1005,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1006,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "CGACTACGATCGACTACGATCGACTACGAT".into(),
        SequenceOperationVariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_different(
        &a.sequence_operation,
        &b.sequence_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(result == false);
}

#[test]
fn test_is_diffable_4() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // INS:chr1:1002
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        1002,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1003,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "CGATC".into(),
        SequenceOperationVariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_different(
        &a.sequence_operation,
        &b.sequence_operation,
        0.05f32,
        0.95f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(result == false);
}

#[test]
fn test_is_diffable_5() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(2, go_1);

    // INS:chr1:1002
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        1002,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1003,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "A".into(),
        SequenceOperationVariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(2, go_2);

    let result: bool = is_different(
        &a.sequence_operation,
        &b.sequence_operation,
        0.05f32,
        0.95f32,
        2000,
        1000,
        1000,
        true
    );

    assert!(result == true);
}

#[test]
fn test_split_variant_records_by_chromosome_1() {
    // INS:chr1:1000
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGACTACGATCGACTACGATCGACT".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(1, go_1);

    // INS:chr2:1002
    let go_2: SequenceOperation = SequenceOperation::new(
        1,
        1002,
        Strand::Forward,
        SequenceOperationType::Downstream,
        1,
        1003,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "CGATC".into(),
        SequenceOperationVariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(1, go_2);

    let mut variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    variant_records.push(Arc::new(a));
    variant_records.push(Arc::new(b));
    let variant_records_map: HashMap<(u16,u16),Vec<Arc<VariantRecord>>> = split_variant_records_by_chromosome(
        variant_records,
        1
    );

    assert!(variant_records_map.get(&(0,0)).unwrap().len() == 1);
    assert!(variant_records_map.get(&(1,1)).unwrap().len() == 1);
}

#[test]
fn test_sweep_clusters_1() {
    // DEL:chr1:1001-1100
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1000,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(1, go_1);

    // DEL:chr1:990-1150
    let go_2: SequenceOperation = SequenceOperation::new(
        0,
        989,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1151,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "".into(),
        SequenceOperationVariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(1, go_2);

    // INS:chr1:1200-1200
    let go_3: SequenceOperation = SequenceOperation::new(
        0,
        1200,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGTAGCTGACGTACATATACTGACC".into(),
        SequenceOperationVariantType::Insertion
    );
    let c: VariantRecord = VariantRecord::new(1, go_3);

    // SNV:chr1:1300
    let go_4: SequenceOperation = SequenceOperation::new(
        0,
        1299,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1301,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "T".into(),
        SequenceOperationVariantType::SingleNucleotideVariant
    );
    let d: VariantRecord = VariantRecord::new(1, go_4);

    let variant_records: Vec<Arc<VariantRecord>> = vec![Arc::new(a), Arc::new(b), Arc::new(c), Arc::new(d)];
    let variant_record_clusters: Vec<VariantRecordCluster> = sweep_clusters(
        variant_records,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        2
    );

    assert!(variant_record_clusters.len() == 3);
}

#[test]
fn test_sweep_clusters_2() {
    // INS:chr1:1200-1200
    let go_1: SequenceOperation = SequenceOperation::new(
        0,
        1200,
        Strand::Forward,
        SequenceOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        SequenceOperationType::Upstream,
        "ACGATCGTAGCTGACGTACATATACTGACC".into(),
        SequenceOperationVariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(1, go_1);

    let variant_records: Vec<Arc<VariantRecord>> = vec![Arc::new(a)];
    let variant_record_clusters: Vec<VariantRecordCluster> = sweep_clusters(
        variant_records,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        2
    );

    assert!(variant_record_clusters.len() == 1);
}
