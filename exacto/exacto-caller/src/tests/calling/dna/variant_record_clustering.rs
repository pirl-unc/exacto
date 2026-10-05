use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use super::*;


/// One read of scga-mini-dna-001, aligned forward, holds the SNV chr17:7,674,225 C>A at base
/// quality 71 and no other event within 100 bases of it (its nearest, a sequencing-error
/// deletion, is at 7,674,361). Its records there are the one SNV, and they cluster into one call
/// of that read.
#[test]
fn cluster_dna_variant_records_returns_one_variant_call_for_scga_mini_dna_001() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    let read_name: &str = "scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/7/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let record: &bam::Record = records_map.get(&read_id).unwrap().first().unwrap();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &vec![Arc::new(record.clone())]
    );

    let alignment: AlignmentModel = alignment.clone();

    // The read's records within 100 bases of the site.
    let chromosome_id: u16 = *chromosome_names_map.get_by_left("chr17").unwrap();
    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(30, 30, 0).call(&alignment)
        .into_iter()
        .filter(|variant_record| variant_record.get_chromosome_1() == chromosome_id
            && (7_674_125..=7_674_325).contains(&variant_record.get_position_1()))
        .collect();

    let variant_records_rc: Vec<Arc<VariantRecord>> = variant_records
        .iter()
        .cloned()
        .map(Arc::new)
        .collect();

    assert!(variant_records.len() == 1);
    assert_eq!(&*variant_records[0].get_graph_operation().as_boxed_str(), "0:7674224:+:D:0:7674226:+:U:A:1:SNV");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
        variant_records_rc,
        0..=u32::MAX,
        &fasta_map,
        &chromosome_names_map,
        1,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        options.calling.bkpt_rescue,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion,
        options.calling.poa_match_score,
        options.calling.poa_mismatch_score,
        options.calling.poa_gap_open_score,
        options.calling.poa_gap_extend_score
    );

    assert!(variant_calls.len() == 1);
    assert_eq!(&*variant_calls[0].get_consensus_graph_operation().as_boxed_str(), "0:7674224:+:D:0:7674226:+:U:A:1:SNV");
    assert_eq!(variant_calls[0].get_read_ids(), vec![read_id]);
}

/// One read of scga-mini-dna-002, aligned forward, holds the 12-base insertion CAGGCGGATGGG after
/// chr17:7,674,224, every inserted base at quality 30 or more, and no other event within 100
/// bases of it (its nearest, a sequencing-error deletion, is at 7,674,361). Its records there
/// are the one insertion, and they cluster into one call of that read. The read's records
/// farther away are pbsim3 sequencing errors and are left out.
#[test]
fn cluster_dna_variant_records_returns_one_variant_call_for_scga_mini_dna_002() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    let read_name: &str = "scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/15/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let record: &bam::Record = records_map.get(&read_id).unwrap().first().unwrap();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &vec![Arc::new(record.clone())]
    );

    let alignment: AlignmentModel = alignment.clone();

    // The read's records within 100 bases of the site.
    let chromosome_id: u16 = *chromosome_names_map.get_by_left("chr17").unwrap();
    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(30, 30, 0).call(&alignment)
        .into_iter()
        .filter(|variant_record| variant_record.get_chromosome_1() == chromosome_id
            && (7_674_125..=7_674_325).contains(&variant_record.get_position_1()))
        .collect();

    let variant_records_rc: Vec<Arc<VariantRecord>> = variant_records
        .iter()
        .cloned()
        .map(Arc::new)
        .collect();

    assert!(variant_records.len() == 1);
    assert_eq!(&*variant_records[0].get_graph_operation().as_boxed_str(), "0:7674224:+:D:0:7674225:+:U:CAGGCGGATGGG:12:INS");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
        variant_records_rc,
        0..=u32::MAX,
        &fasta_map,
        &chromosome_names_map,
        1,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        options.calling.bkpt_rescue,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion,
        options.calling.poa_match_score,
        options.calling.poa_mismatch_score,
        options.calling.poa_gap_open_score,
        options.calling.poa_gap_extend_score
    );

    assert!(variant_calls.len() == 1);
    assert_eq!(&*variant_calls[0].get_consensus_graph_operation().as_boxed_str(), "0:7674224:+:D:0:7674225:+:U:CAGGCGGATGGG:12:INS");
    assert_eq!(variant_calls[0].get_read_ids(), vec![read_id]);
}

#[test]
fn meet_dna_clustering_criteria_returns_false_for_different_variant_types() {
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
        "ACGATCGACT".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_1
    );

    // DEL:chr1:1001-1100
    let go_2: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        2,
        3,
        go_2
    );

    let result: bool = meet_dna_clustering_criteria(
        &a,
        &b,
        0.5f64,
        0.5f64,
        1000,
        0.01f64
    );

    assert!(result == false);
}

#[test]
fn meet_dna_clustering_criteria_returns_true_for_shifted_similar_insertions() {
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
        "ACGATCGACT".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_1
    );

    // INS:chr1:1001
    let go_2: GraphOperation = GraphOperation::new(
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1002,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CGATCGACTC".into(),
        VariantType::Insertion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        2,
        3,
        go_2
    );

    let result: bool = meet_dna_clustering_criteria(
        &a,
        &b,
        0.5f64,
        0.5f64,
        1000,
        0.01f64
    );

    assert!(result == true);
}

#[test]
fn meet_dna_clustering_criteria_returns_true_for_proximal_deletions_of_similar_size() {
    // DEL:chr1:1001-1100
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        0,
        0,
        go_1
    );

    // DEL:chr1:999-1110
    let go_2: GraphOperation = GraphOperation::new(
        0,
        998,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1111,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let b: VariantRecord = VariantRecord::new(
        2,
        10,
        10,
        go_2
    );

    let result: bool = meet_dna_clustering_criteria(
        &a,
        &b,
        0.5f64,
        0.5f64,
        1000,
        0.01f64
    );

    assert!(result == true);
}

#[test]
fn meet_dna_clustering_criteria_returns_true_for_proximal_interchromosomal_breakpoints() {
    // TRA:chr1:1001-chr2:2001
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Downstream,
        1,
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
        200,
        200,
        go_2
    );

    let result: bool = meet_dna_clustering_criteria(
        &a,
        &b,
        0.5f64,
        0.5f64,
        1000,
        0.01f64
    );

    assert!(result == true);
}

#[test]
fn split_variant_records_separates_records_on_different_chromosomes() {
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

    // INS:chr2:1002
    let go_2: GraphOperation = GraphOperation::new(
        1,
        1002,
        Strand::Forward,
        GraphOperationType::Downstream,
        1,
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

    let mut variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    variant_records.push(Arc::new(a));
    variant_records.push(Arc::new(b));
    let variant_records_map: HashMap<(u16, u16, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> = split_variant_records(
        variant_records,
        true,
        1
    );

    assert!(variant_records_map.get(&(0,0,VariantType::Insertion, GraphOperationType::Downstream, GraphOperationType::Upstream)).unwrap().len() == 1);
    assert!(variant_records_map.get(&(1,1,VariantType::Insertion, GraphOperationType::Downstream, GraphOperationType::Upstream)).unwrap().len() == 1);
}

#[test]
fn cluster_non_breakpoint_dna_variant_records_merges_only_the_proximal_deletions() {
    // DEL:chr1:1001-1100
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1101,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        5,
        5,
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
        10,
        10,
        go_2
    );

    // INS:chr1:1200-1200
    let go_3: GraphOperation = GraphOperation::new(
        0,
        1200,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGTAGCTGACGTACATATACTGACC".into(),
        VariantType::Insertion
    );
    let c: VariantRecord = VariantRecord::new(
        1,
        31,
        60,
        go_3
    );

    // SNV:chr1:1300
    let go_4: GraphOperation = GraphOperation::new(
        0,
        1299,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1301,
        Strand::Forward,
        GraphOperationType::Upstream,
        "T".into(),
        VariantType::SingleNucleotideVariant
    );
    let d: VariantRecord = VariantRecord::new(
        1,
        5,
        5,
        go_4
    );

    let variant_records: Vec<Arc<VariantRecord>> = vec![Arc::new(a), Arc::new(b), Arc::new(c), Arc::new(d)];

    let variant_record_clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_dna_variant_records(
        &variant_records,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        2
    );

    assert!(variant_record_clusters.len() == 3);
}

#[test]
fn cluster_non_breakpoint_dna_variant_records_returns_one_cluster_for_single_record() {
    // INS:chr1:1200-1200
    let go_1: GraphOperation = GraphOperation::new(
        0,
        1200,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGTAGCTGACGTACATATACTGACC".into(),
        VariantType::Insertion
    );
    let a: VariantRecord = VariantRecord::new(
        1,
        31,
        60,
        go_1
    );

    let variant_records: Vec<Arc<VariantRecord>> = vec![Arc::new(a)];

    let variant_record_clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_dna_variant_records(
        &variant_records,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        2
    );

    assert!(variant_record_clusters.len() == 1);
}

#[test]
fn cluster_non_breakpoint_dna_variant_records_is_deterministic_across_repeated_calls() {
    let record = |chromosome: u16, position: u32, sequence: &str| -> Arc<VariantRecord> {
        let graph_operation: GraphOperation = GraphOperation::new(
            chromosome,
            position,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome,
            position + 1,
            Strand::Forward,
            GraphOperationType::Upstream,
            sequence.into(),
            VariantType::Insertion
        );
        Arc::new(VariantRecord::new(0, 6, 10, graph_operation))
    };

    let variant_records: Vec<Arc<VariantRecord>> = (0..24u32)
        .map(|i| record(0, 1000 + i * 10, "CGATC"))
        .collect();

    let expected_map = split_variant_records(
        variant_records.clone(),
        true,
        1
    );

    let expected_clusters = cluster_non_breakpoint_dna_variant_records(
        &variant_records,
        0.5,
        0.5,
        1000,
        0.01,
        1
    );

    let num_jobs: usize = 100_000;
    let mut outcomes: Vec<(usize, usize)> = Vec::with_capacity(num_jobs);
    for _ in 0..num_jobs {
        let map = split_variant_records(
            variant_records.clone(),
            true,
            1
        );
        let clusters = cluster_non_breakpoint_dna_variant_records(
            &variant_records,
            0.5,
            0.5,
            1000,
            0.01,
            1
        );
        outcomes.push((map.len(), clusters.len()));
    }

    assert_eq!(outcomes.len(), 100_000);

    for (num_keys, num_clusters) in outcomes.iter() {
        assert_eq!(*num_keys, expected_map.len(), "split_variant_records differed inside a pool");
        assert_eq!(*num_clusters, expected_clusters.len(), "sweep_clusters differed inside a pool");
    }
}

/// One unresolved clip: a `Breakpoint` whose second operation is `Noop`.
fn softclip_breakend(read_id: usize, position: u32) -> Arc<VariantRecord> {
    Arc::new(VariantRecord::new(
        read_id,
        0,
        10,
        GraphOperation::new(
            0,
            position,
            Strand::Forward,
            GraphOperationType::Downstream,
            0,
            position,
            Strand::Forward,
            GraphOperationType::Noop,
            "ACGTACGTAC".to_string().into_boxed_str(),
            VariantType::Breakpoint
        )
    ))
}

#[test]
fn cluster_breakpoint_dna_variant_records_pools_nearby_orphan_clips_as_one_breakpoint() {
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(
        vec![softclip_breakend(1, 1000), softclip_breakend(2, 1050)],
        100,
        0.5f64,
        0.5f64,
        0.01f64,
        1
    );
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].get_variant_records().len(), 2);

    // The rebuilt insertions sit 50 bp apart, outside the size-scaled tolerance,
    // so the orphan pool is kept as a breakpoint cluster rather than collapsed
    // to a terminal insertion.
    assert!(clusters[0]
        .get_variant_records()
        .iter()
        .all(|vr| vr.get_variant_type() == &VariantType::Breakpoint));
}

#[test]
fn cluster_breakpoint_dna_variant_records_rebuilds_agreeing_orphan_clips_as_one_insertion() {
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(
        vec![softclip_breakend(1, 1000), softclip_breakend(2, 1000)],
        100,
        0.5f64,
        0.5f64,
        0.001f64,
        1
    );
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].get_variant_records().len(), 2);

    // Both tails carry the same inserted allele at the same anchor, so the pool
    // is a terminal insertion, not a junction.
    assert!(clusters[0]
        .get_variant_records()
        .iter()
        .all(|vr| vr.get_variant_type() == &VariantType::Insertion));
}

#[test]
fn cluster_breakpoint_dna_variant_records_absorbs_orphan_clips_into_resolved_split_read() {
    // A resolved junction anchors both tails: 1000 and 1050 each sit within the
    // tolerance of the breakend at 1020, so all three land in one pool and the
    // clips become read support for the junction.
    let split_read: Arc<VariantRecord> = Arc::new(VariantRecord::new(
        3,
        0,
        10,
        GraphOperation::new(
            0,
            1020,
            Strand::Forward,
            GraphOperationType::Downstream,
            0,
            5000,
            Strand::Reverse,
            GraphOperationType::Downstream,
            "".to_string().into_boxed_str(),
            VariantType::Breakpoint
        )
    ));
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(
        vec![softclip_breakend(1, 1000), softclip_breakend(2, 1050), split_read],
        100,
        0.5f64,
        0.5f64,
        0.001f64,
        1
    );
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].get_variant_records().len(), 3);
    assert!(clusters[0]
        .get_variant_records()
        .iter()
        .any(|vr| vr.get_read_id() == 3));
}

/// A translocation chr0:1000 -> chr1:5000 held by three split reads, and a fourth read clipped
/// at the chr0 breakend. The clip's own chromosome 2 is chr0 (it has no mate), yet it is read
/// support for the translocation: the cluster takes its chromosome 2 from the split reads, and
/// the call holds all four reads.
#[test]
fn cluster_breakpoint_dna_variant_records_keeps_a_clip_in_a_translocation_cluster() {
    let mut records: Vec<Arc<VariantRecord>> = Vec::new();
    for read_id in 1..=3 {
        records.push(Arc::new(VariantRecord::new(
            read_id,
            500,
            501,
            GraphOperation::new(
                0, 1_000, Strand::Forward, GraphOperationType::Downstream,
                1, 5_000, Strand::Forward, GraphOperationType::Upstream,
                "".into(),
                VariantType::Translocation
            )
        )));
    }
    records.push(Arc::new(VariantRecord::new(
        4,
        500,
        539,
        GraphOperation::new(
            0, 1_000, Strand::Forward, GraphOperationType::Downstream,
            0, 1_000, Strand::Forward, GraphOperationType::Noop,
            "ACGTTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAGTT".into(),
            VariantType::Breakpoint
        )
    )));

    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(records, 1_000, 0.5, 0.5, 0.01, 1);

    assert_eq!(clusters.len(), 1);
    let variant_records: HashSet<VariantRecord> = clusters[0].get_variant_records().iter().map(|record| (**record).clone()).collect();
    let variant_call: VariantCall = VariantCall::from_variant_records(0, variant_records, 0, 4, 6, 2);
    assert_eq!(variant_call.get_read_ids(), vec![1, 2, 3, 4]);
    let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
    assert_eq!(
        (operation.get_chromosome_1(), operation.get_position_1(), operation.get_chromosome_2(), operation.get_position_2()),
        (0, 1_000, 1, 5_000)
    );
    assert_eq!(*operation.get_variant_type(), VariantType::Translocation);
}

/// A junction 1000 -> 5000 held by three split reads, two reads clipped at its 1000 breakend
/// (Downstream) and two at its 5000 breakend (Upstream). A clip at 5000 has no resolved record in
/// its own position 1 group, but its breakend is the junction's side 2, so all seven reads
/// support one junction. A clip at 5000 that points the other way (Downstream) is not this
/// junction's breakend and stays out, and so does the 6-base clip that ends a read 306 bases
/// past it: a breakend so short a clip can only lie a few bases off.
#[test]
fn cluster_breakpoint_dna_variant_records_joins_clips_at_either_breakend_of_a_junction() {
    let mut records: Vec<Arc<VariantRecord>> = Vec::new();
    for (read_id, position_1, position_2) in [(1, 1_000, 5_000), (2, 1_000, 5_000), (3, 1_001, 5_001)] {
        records.push(Arc::new(VariantRecord::new(
            read_id,
            500,
            501,
            GraphOperation::new(
                0, position_1, Strand::Forward, GraphOperationType::Downstream,
                0, position_2, Strand::Forward, GraphOperationType::Upstream,
                "".into(),
                VariantType::Breakpoint
            )
        )));
    }
    for (read_id, position, operation, sequence) in [
        (4, 1_000, GraphOperationType::Downstream, "ACGTTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAG"),
        (5, 1_000, GraphOperationType::Downstream, "ACGTTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAGTTCAGCATGGACTTGCACGAAT"),
        (6, 5_000, GraphOperationType::Upstream, "GGCATCCGATAAGTCGTACCTGAG"),
        (7, 5_000, GraphOperationType::Upstream, "TTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAGTTCAGCATGGACTTG"),
        (8, 5_000, GraphOperationType::Downstream, "CCGATAAGTCGTACCTGAGTTCAGCATGGA"),
        (9, 5_306, GraphOperationType::Upstream, "ACGTTG")
    ] {
        records.push(Arc::new(VariantRecord::new(
            read_id,
            0,
            sequence.len() as u32 - 1,
            GraphOperation::new(
                0, position, Strand::Forward, operation,
                0, position, Strand::Forward, GraphOperationType::Noop,
                sequence.into(),
                VariantType::Breakpoint
            )
        )));
    }

    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(records, 1_000, 0.5, 0.5, 0.01, 1);

    let mut read_ids: Vec<Vec<usize>> = clusters
        .iter()
        .map(|cluster| {
            let mut read_ids: Vec<usize> = cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect();
            read_ids.sort();
            read_ids
        })
        .collect();
    read_ids.sort();
    assert_eq!(read_ids, vec![vec![1, 2, 3, 4, 5, 6, 7], vec![8], vec![9]]);

    let junction: &VariantRecordCluster = clusters.iter().find(|cluster| cluster.get_variant_records().len() == 7).unwrap();
    let variant_records: HashSet<VariantRecord> = junction.get_variant_records().iter().map(|record| (**record).clone()).collect();
    let variant_call: VariantCall = VariantCall::from_variant_records(0, variant_records, 0, 4, 6, 2);
    let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
    assert_eq!((operation.get_position_1(), operation.get_position_2()), (1_000, 5_000));
    assert_eq!(variant_call.get_num_reads(), 7);
}

/// The same for a translocation: a read clipped at the chr1 breakend of chr0:1000 -> chr1:5000
/// has chromosome 1 = chr1, the junction's chromosome 2, and joins it.
#[test]
fn cluster_breakpoint_dna_variant_records_joins_a_clip_at_the_partner_breakend_of_a_translocation() {
    let mut records: Vec<Arc<VariantRecord>> = Vec::new();
    for read_id in 1..=3 {
        records.push(Arc::new(VariantRecord::new(
            read_id,
            500,
            501,
            GraphOperation::new(
                0, 1_000, Strand::Forward, GraphOperationType::Downstream,
                1, 5_000, Strand::Forward, GraphOperationType::Upstream,
                "".into(),
                VariantType::Translocation
            )
        )));
    }
    records.push(Arc::new(VariantRecord::new(
        4,
        0,
        29,
        GraphOperation::new(
            1, 5_002, Strand::Forward, GraphOperationType::Upstream,
            1, 5_002, Strand::Forward, GraphOperationType::Noop,
            "CCGATAAGTCGTACCTGAGTTCAGCATGGA".into(),
            VariantType::Breakpoint
        )
    )));

    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(records, 1_000, 0.5, 0.5, 0.01, 1);

    assert_eq!(clusters.len(), 1);
    let variant_records: HashSet<VariantRecord> = clusters[0].get_variant_records().iter().map(|record| (**record).clone()).collect();
    let variant_call: VariantCall = VariantCall::from_variant_records(0, variant_records, 0, 4, 6, 2);
    assert_eq!(variant_call.get_read_ids(), vec![1, 2, 3, 4]);
    assert_eq!(*variant_call.get_consensus_graph_operation().get_variant_type(), VariantType::Translocation);
}

/// scga-mini-dna-001 covers chr17:7,668,421-7,687,490 and nothing else (chr18 carries no reads).
/// The depths, as `samtools depth -a` counts them: 0 at 7,668,420, 10 at 7,668,421, 66 at
/// 7,674,224 and 7,674,226, and 65 at 7,674,225, where one read holds a deletion. A local
/// variant's denominator is the depth at the base(s) it replaces: one base for an SNV, the
/// deepest of the replaced run for an MNV, the deepest of the two flanks for an insertion.
#[test]
fn get_dna_variant_position_total_depth_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let positions: HashMap<ReferenceChromosomeName, Vec<ReferencePosition>> = HashMap::from([
        ("chr17".into(), (7_660_000..=7_687_500).collect()),
        ("chr18".into(), vec![5_170_100])
    ]);
    let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &format!("{bam_file}.bai"), &positions, 1_000);

    // SNV at 7,674,225, inside the covered span: its own base, not its deeper flanks.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_674_224, &GraphOperationType::Downstream,
            "chr17", 7_674_226, &GraphOperationType::Upstream,
            "A"
        ),
        65
    );

    // SNV at 7,660,001, outside it.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_660_000, &GraphOperationType::Downstream,
            "chr17", 7_660_002, &GraphOperationType::Upstream,
            "A"
        ),
        0
    );

    // MNV replacing 7,668,420-7,668,421: the first base is uncovered, the second is not.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_419, &GraphOperationType::Downstream,
            "chr17", 7_668_422, &GraphOperationType::Upstream,
            "AA"
        ),
        10
    );

    // Insertion between 7,668,420 and 7,668,421: one flank uncovered, one covered.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_420, &GraphOperationType::Downstream,
            "chr17", 7_668_421, &GraphOperationType::Upstream,
            "A"
        ),
        10
    );
}

/// A deletion's denominator is the deeper of its two flanking bases, whichever side is
/// deeper, and never the bases it removes. In scga-mini-dna-001 7,668,415 and 7,687,495 are
/// uncovered, and 7,668,425 and 7,687,485 are both at depth 12 (`samtools depth -a`).
#[test]
fn get_dna_variant_position_total_depth_takes_the_deeper_end_of_a_deletion() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let positions: HashMap<ReferenceChromosomeName, Vec<ReferencePosition>> = HashMap::from([
        ("chr17".into(), (7_660_000..=7_687_500).collect()),
        ("chr18".into(), vec![5_170_100])
    ]);
    let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &format!("{bam_file}.bai"), &positions, 1_000);

    // Shallow end first: 7,668,415 is uncovered, 7,668,425 is covered.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_415, &GraphOperationType::Downstream,
            "chr17", 7_668_425, &GraphOperationType::Upstream,
            ""
        ),
        12
    );

    // Deep end first: 7,687,485 is covered, 7,687,495 is not.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_687_485, &GraphOperationType::Downstream,
            "chr17", 7_687_495, &GraphOperationType::Upstream,
            ""
        ),
        12
    );

    // Both ends uncovered with the whole covered span deleted between them: the interior is
    // never consulted.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_415, &GraphOperationType::Downstream,
            "chr17", 7_687_495, &GraphOperationType::Upstream,
            ""
        ),
        0
    );
}

/// A breakend's denominator is the deeper of its two sides, on one chromosome or two. In
/// scga-mini-dna-001 7,668,415 is uncovered, 7,668,425 is at depth 12 and 7,674,224 at 66
/// (`samtools depth -a`); chr18 carries no reads.
#[test]
fn get_dna_variant_position_total_depth_takes_the_deeper_end_of_a_breakend() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let positions: HashMap<ReferenceChromosomeName, Vec<ReferencePosition>> = HashMap::from([
        ("chr17".into(), (7_660_000..=7_687_500).collect()),
        ("chr18".into(), vec![5_170_100])
    ]);
    let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &format!("{bam_file}.bai"), &positions, 1_000);

    // Same chromosome, both sides Downstream (an inversion-shaped junction), straddling the
    // coverage edge.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_415, &GraphOperationType::Downstream,
            "chr17", 7_668_425, &GraphOperationType::Downstream,
            ""
        ),
        12
    );

    // Across chromosomes: chr18 carries no reads, so the chr17 side decides either way round.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_674_224, &GraphOperationType::Downstream,
            "chr18", 5_170_100, &GraphOperationType::Upstream,
            ""
        ),
        66
    );
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr18", 5_170_100, &GraphOperationType::Upstream,
            "chr17", 7_674_224, &GraphOperationType::Downstream,
            ""
        ),
        66
    );

    // Neither side covered.
    assert_eq!(
        get_dna_variant_position_total_depth(
            &read_depths,
            "chr17", 7_668_415, &GraphOperationType::Downstream,
            "chr18", 5_170_100, &GraphOperationType::Upstream,
            ""
        ),
        0
    );
}


/// Three reads hold A, two hold T and one holds C at chr17:7,674,225. Each base is a call of its
/// own: a read that holds another base is not support for the consensus base.
#[test]
fn cluster_dna_variant_records_keeps_snv_bases_at_one_position_apart() {
    let fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    let variant_records: Vec<Arc<VariantRecord>> = [(1, "A"), (2, "A"), (3, "A"), (4, "T"), (5, "T"), (6, "C")]
        .iter()
        .map(|&(read_id, base)| {
            Arc::new(VariantRecord::new(
                read_id,
                100,
                100,
                GraphOperation::new(
                    0,
                    7_674_224,
                    Strand::Forward,
                    GraphOperationType::Downstream,
                    0,
                    7_674_226,
                    Strand::Forward,
                    GraphOperationType::Upstream,
                    base.into(),
                    VariantType::SingleNucleotideVariant
                )
            ))
        })
        .collect();
    let options: DNAVariantCallingOptions = DNAVariantCallingOptions::DEFAULT;

    let variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
        variant_records,
        0..=u32::MAX,
        &fasta_map,
        &chromosome_names_map,
        1,
        options.min_size_proportion,
        options.max_ins_norm_edit_distance,
        options.max_clustering_distance,
        0.01,
        options.bkpt_rescue,
        options.bkpt_rescue_min_ins_len,
        options.bkpt_rescue_max_ins_len,
        options.bkpt_rescue_search_distance,
        options.bkpt_rescue_realignment_gap_open_score,
        options.bkpt_rescue_realignment_gap_extend_score,
        options.bkpt_rescue_realignment_k,
        options.bkpt_rescue_realignment_band_width,
        options.bkpt_rescue_realignment_min_score_fraction,
        options.bkpt_rescue_realignment_min_query_coverage,
        options.bkpt_rescue_realignment_min_span_proportion,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    let mut alleles: Vec<(String, Vec<usize>)> = variant_calls
        .iter()
        .map(|variant_call| (variant_call.get_consensus_graph_operation().get_standardized_sequence(), variant_call.get_read_ids()))
        .collect();
    alleles.sort();
    assert_eq!(
        alleles,
        vec![("A".to_string(), vec![1, 2, 3]), ("C".to_string(), vec![6]), ("T".to_string(), vec![4, 5])]
    );
}


/// Two reads hold the inverted copy of chr17:7,673,001-7,673,500 as an insertion after 7,673,500,
/// and a third splits there into the copy (the fold-back at 7,673,500). Breakpoint rescue places
/// the insertions on the reference and retypes each as the fold-back and the far junction
/// 7,673,001 -> 7,673,501, whose position 1 lies 500 bases before the insertions. A call starts
/// at the records it was built from, the insertions for a retyped one, so the window whose chunk
/// holds 7,673,500 makes both calls, even when 7,673,001 lies in the chunk before it.
#[test]
fn cluster_dna_variant_records_starts_a_retyped_insertion_at_the_insertion() {
    let fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    let inverted_copy: String = reverse_complement(&fasta_map.get_sequence("chr17", 7_673_001, 7_673_500).to_ascii_uppercase()).to_string();
    let mut variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    for read_id in [1, 2] {
        variant_records.push(Arc::new(VariantRecord::new(
            read_id,
            1_000,
            1_499,
            GraphOperation::new(
                0, 7_673_500, Strand::Forward, GraphOperationType::Downstream,
                0, 7_673_501, Strand::Forward, GraphOperationType::Upstream,
                inverted_copy.clone().into(),
                VariantType::Insertion
            )
        )));
    }
    variant_records.push(Arc::new(VariantRecord::new(
        3,
        1_000,
        1_001,
        GraphOperation::new(
            0, 7_673_500, Strand::Forward, GraphOperationType::Downstream,
            0, 7_673_500, Strand::Reverse, GraphOperationType::Downstream,
            "".into(),
            VariantType::Breakpoint
        )
    )));
    let options: DNAVariantCallingOptions = DNAVariantCallingOptions::DEFAULT;

    let calls = |owned_positions: std::ops::RangeInclusive<u32>| -> Vec<(Box<str>, Vec<usize>)> {
        let mut calls: Vec<(Box<str>, Vec<usize>)> = cluster_dna_variant_records(
            variant_records.clone(),
            owned_positions,
            &fasta_map,
            &chromosome_names_map,
            1,
            options.min_size_proportion,
            options.max_ins_norm_edit_distance,
            options.max_clustering_distance,
            0.01,
            true,
            options.bkpt_rescue_min_ins_len,
            options.bkpt_rescue_max_ins_len,
            options.bkpt_rescue_search_distance,
            options.bkpt_rescue_realignment_gap_open_score,
            options.bkpt_rescue_realignment_gap_extend_score,
            options.bkpt_rescue_realignment_k,
            options.bkpt_rescue_realignment_band_width,
            options.bkpt_rescue_realignment_min_score_fraction,
            options.bkpt_rescue_realignment_min_query_coverage,
            options.bkpt_rescue_realignment_min_span_proportion,
            options.poa_match_score,
            options.poa_mismatch_score,
            options.poa_gap_open_score,
            options.poa_gap_extend_score
        )
            .iter()
            .map(|variant_call| (variant_call.get_consensus_graph_operation().as_boxed_str(), variant_call.get_read_ids()))
            .collect();
        calls.sort();
        calls
    };

    let all_calls: Vec<(Box<str>, Vec<usize>)> = calls(0..=u32::MAX);
    assert_eq!(all_calls, vec![
        ("0:7673001:-:U:0:7673501:+:U::0:BND".into(), vec![1, 2]),
        ("0:7673500:+:D:0:7673500:-:D::0:BND".into(), vec![1, 2, 3])
    ]);
    assert_eq!(calls(7_673_101..=7_674_000), all_calls);
    assert!(calls(0..=7_673_499).is_empty());
    assert!(calls(7_673_501..=u32::MAX).is_empty());
}


/// A clip on the first base of a contig is an insertion at (0, 1), one on the last base of a
/// contig of length L an insertion at (L, L + 1). The flank off the contig is read at the base
/// beside it.
#[test]
fn get_dna_variant_position_total_depth_reads_a_flank_off_the_contig_at_the_base_beside_it() {
    // Contig lengths from the header (chr18: 10,000,000 bases, no reads); the counts are set by hand.
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let mut read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &format!("{bam_file}.bai"), &HashMap::new(), 1_000);
    read_depths.insert("chr18", 1, 11, 0, 0);
    read_depths.insert("chr18", 2, 30, 0, 0);
    read_depths.insert("chr18", 10_000_000, 13, 0, 0);

    // Insertion before the first base and after the last base.
    let first: u32 = get_dna_variant_position_total_depth(
        &read_depths, "chr18", 0, &GraphOperationType::Downstream, "chr18", 1, &GraphOperationType::Upstream, "ACGT"
    );
    let last: u32 = get_dna_variant_position_total_depth(
        &read_depths, "chr18", 10_000_000, &GraphOperationType::Downstream, "chr18", 10_000_001, &GraphOperationType::Upstream, "ACGT"
    );
    assert_eq!((first, last), (11, 13));

    // Deletion of the first base, SNV on the last base: the flank inside the contig is read.
    let deletion: u32 = get_dna_variant_position_total_depth(
        &read_depths, "chr18", 0, &GraphOperationType::Downstream, "chr18", 2, &GraphOperationType::Upstream, ""
    );
    let snv: u32 = get_dna_variant_position_total_depth(
        &read_depths, "chr18", 9_999_999, &GraphOperationType::Downstream, "chr18", 10_000_001, &GraphOperationType::Upstream, "A"
    );
    assert_eq!((deletion, snv), (30, 13));
}


/// Read 1 holds two insertions of 50 bases 30 bases apart (one expansion the aligner wrote as
/// two); reads 2 and 3 hold one between them. The call has three reads, not four.
#[test]
fn cluster_dna_variant_records_counts_a_read_with_two_records_in_a_call_once() {
    let fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    let sequence: String = "ACGTTGCA".repeat(7)[..50].to_string();
    let variant_records: Vec<Arc<VariantRecord>> = [(1, 1000, 7_674_200), (1, 1080, 7_674_230), (2, 2000, 7_674_215), (3, 3000, 7_674_215)]
        .iter()
        .map(|&(read_id, read_position, anchor)| {
            Arc::new(VariantRecord::new(
                read_id,
                read_position,
                read_position + 50,
                GraphOperation::new(
                    0,
                    anchor,
                    Strand::Forward,
                    GraphOperationType::Downstream,
                    0,
                    anchor + 1,
                    Strand::Forward,
                    GraphOperationType::Upstream,
                    sequence.as_str().into(),
                    VariantType::Insertion
                )
            ))
        })
        .collect();
    let options: DNAVariantCallingOptions = DNAVariantCallingOptions::DEFAULT;

    let mut variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
        variant_records,
        0..=u32::MAX,
        &fasta_map,
        &chromosome_names_map,
        1,
        options.min_size_proportion,
        options.max_ins_norm_edit_distance,
        options.max_clustering_distance,
        0.01,
        options.bkpt_rescue,
        options.bkpt_rescue_min_ins_len,
        options.bkpt_rescue_max_ins_len,
        options.bkpt_rescue_search_distance,
        options.bkpt_rescue_realignment_gap_open_score,
        options.bkpt_rescue_realignment_gap_extend_score,
        options.bkpt_rescue_realignment_k,
        options.bkpt_rescue_realignment_band_width,
        options.bkpt_rescue_realignment_min_score_fraction,
        options.bkpt_rescue_realignment_min_query_coverage,
        options.bkpt_rescue_realignment_min_span_proportion,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(variant_calls.len(), 1);
    assert_eq!(variant_calls[0].get_variant_records().len(), 4);
    assert_eq!(variant_calls[0].get_read_ids(), vec![1, 2, 3]);
    // The pipeline sets the total depth from the reads after clustering; 30 reads cover the site.
    variant_calls[0].set_total_depth(30);
    assert_eq!(variant_calls[0].get_alternate_allele_fraction(), 3.0 / 30.0);
}
