use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;
use bimap::BiMap;

use crate::prelude::*;


#[test]
fn test_compute_min_read_support_1() {
    let (min_read_support, f1, recall, fpr, precision) = compute_min_read_support(
        28,
        0.5f64,
        1e-3,
        1e-2,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support, 6);

    let (min_read_support, f1, recall, fpr, precision) = compute_min_read_support(
        30,
        0.5f64,
        1e-3,
        1e-2,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support, 6);

    let (min_read_support, f1, recall, fpr, precision) = compute_min_read_support(
        60,
        0.25f64,
        1e-6,
        1e-2,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support, 9);
}

#[test]
fn test_compute_min_read_support_index_1() {
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        60u64,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index[0][29], 6);
}

#[test]
fn test_compute_min_read_support_index_2() {
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        60u64,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index[4][29], 11);
}

#[test]
fn test_compute_min_read_support_index_3() {
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        60u64,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index[5][29], 13);
}

#[test]
fn test_compute_min_read_support_index_4() {
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        60u64,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        0.99f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index[6][29], 14);
}

#[test]
fn test_get_variant_position_total_depth() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let depths_map = get_bam_depths_map(bam_file, 2);
    let total_depth: u32 = get_variant_position_total_depth(
        &depths_map,
        "chr17",
        7_674_224,
        &GraphOperationType::Downstream,
        "chr17",
        7_674_226,
        &GraphOperationType::Upstream,
        "A"
    );

    assert!(total_depth == 6);
}

#[test]
fn test_has_strand_bias() {
    let strand_bias_exists: bool = has_strand_bias(
        0,
        8,
        12,
        13,
        0.05
    );
    assert_eq!(strand_bias_exists, true);
}

#[test]
fn test_is_repeat_indel_1() {
    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_837,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_841,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);
}

#[test]
fn test_is_repeat_indel_2() {
    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_837,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Upstream,
        "AAAA".into(),
        VariantType::Insertion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        104,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 9);
}

#[test]
fn test_is_repeat_indel_3() {
    let go: GraphOperation = GraphOperation::new(
        0,
        4_329_639,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_329_640,
        Strand::Forward,
        GraphOperationType::Upstream,
        "AT".into(),
        VariantType::Insertion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 8);
}

#[test]
fn test_is_repeat_indel_4() {
    let go: GraphOperation = GraphOperation::new(
        0,
        4_329_640,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_329_641,
        Strand::Forward,
        GraphOperationType::Upstream,
        "TA".into(),
        VariantType::Insertion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 6);
}



#[test]
fn test_is_repeat_indel_5() {
    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_840,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);
}

#[test]
fn test_is_repeat_indel_6() {
    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_836,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);
}

#[test]
fn test_is_repeat_indel_7() {
    let go: GraphOperation = GraphOperation::new(
        0,
        4_330_358,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_330_361,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let (is_repeat, size) = is_repeat_indel(
        &vr,
        &chromosome_names_map,
        &fasta_map
    );

    assert_eq!(is_repeat, false);
}

#[test]
fn test_variant_calling_1() {
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
    let mut variant_call: VariantCall = VariantCall::new(1);
    variant_call.add_variant_record(variant_record_1);
    variant_call.add_variant_record(variant_record_2);
    assert!(variant_call.get_consensus_record().1.len() == 2);
}

#[test]
fn test_variant_calling_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let depths_map: Arc<HashMap<Box<str>, Vec<u32>>> = Arc::new(get_bam_depths_map(bam_file, 2));
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
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
        "chr17",
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    let read_name: &str = "m64012_325382_158010/1/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let record: &bam::Record = records_map.get(&read_id).unwrap().first().unwrap();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &vec![Arc::new(record.clone())]
    );
    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();
    let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(30, 30, AnalyteType::DNA);
    assert!(variant_records.len() == 1);
    let variant_calls: Vec<VariantCall> = cluster_variant_records(
        variant_records.iter().map(|record| Arc::new(record.clone())).collect(),
        &depths_map,
        &chromosome_names_map,
        1,
        0.5f32,
        0.5f32,
        2000,
        1000,
        10000,
        false
    );
    assert!(variant_calls.len() == 1);
}

#[test]
fn test_variant_calling_3() {
    // INS:chr1:1000
    // DEL:chr1:1100-1200
    // SNV:chr1:1205
    // SNV:chr1:1400
    let mut a: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1001
    // SNV:chr1:1205
    let mut b: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1000
    let go_a1: GraphOperation = GraphOperation::new(
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
    let vr_a1: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_a1
    );

    // DEL:chr1:1101-1200
    let go_a2: GraphOperation = GraphOperation::new(
        0,
        1100,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_a2: VariantRecord = VariantRecord::new(
        1,
        2,
        104,
        go_a2
    );

    // SNV:chr1:1205
    let go_a3: GraphOperation = GraphOperation::new(
        0,
        1204,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_a3: VariantRecord = VariantRecord::new(
        1,
        105,
        107,
        go_a3
    );

    // SNV:chr1:1400
    let go_a4: GraphOperation = GraphOperation::new(
        0,
        1399,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1401,
        Strand::Forward,
        GraphOperationType::Upstream,
        "T".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_a4: VariantRecord = VariantRecord::new(
        1,
        108,
        110,
        go_a4
    );

    a.push(vr_a1);
    a.push(vr_a2);
    a.push(vr_a3);
    a.push(vr_a4);

    // INS:chr1:1001
    let go_b1: GraphOperation = GraphOperation::new(
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1002,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CGATCGACTA".into(),
        VariantType::Insertion
    );
    let vr_b1: VariantRecord = VariantRecord::new(
        2,
        111,
        112,
        go_b1
    );

    // SNV:chr1:1205
    let go_b2: GraphOperation = GraphOperation::new(
        0,
        1204,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_b2: VariantRecord = VariantRecord::new(
        2,
        113,
        115,
        go_b2
    );

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
        true,
        false
    );

    assert!(a_only.len() == 2);

    let mut found_del_1: bool = false;
    let mut found_snv_2: bool = false;

    for variant_record in a_only.iter() {
        if variant_record.get_graph_operation_boxed_str() == "0:1100:+:D:0:1201:+:U::0:DEL".into() {
            found_del_1 = true;
        }
        if variant_record.get_graph_operation_boxed_str() == "0:1399:+:D:0:1401:+:U:T:1:SNV".into() {
            found_snv_2 = true;
        }
    }

    assert!(found_del_1);
    assert!(found_snv_2);
}

#[test]
fn test_variant_calling_4() {
    let bam_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let depths_map: Arc<HashMap<Box<str>, Vec<u32>>> = Arc::new(get_bam_depths_map(bam_file, 2));
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
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
        "chr17",
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    let read_name: &str = "m64012_202369_785869/3/ccs";
    let read_id: usize = read_names_map.get_by_left(&read_name.to_string().into_boxed_str()).unwrap().clone();
    let read_sequence: Box<str> = get_primary_alignment_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_primary_alignment_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );
    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();
    let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(
        30,
        30,
        AnalyteType::DNA
    );
    assert!(variant_records.len() == 1);
    let variant_calls: Vec<VariantCall> = cluster_variant_records(
        variant_records.iter().map(|record| Arc::new(record.clone())).collect(),
        &depths_map,
        &chromosome_names_map,
        1,
        0.5f32,
        0.5f32,
        2000,
        1000,
        10000,
        false
    );
    assert!(variant_calls.len() == 1);
}

#[test]
fn test_variant_calling_5() {
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

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        false
    );

    assert!(result == false);
}

#[test]
fn test_variant_calling_6() {
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

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        false
    );

    assert!(result == true);
}

#[test]
fn test_variant_calling_7() {
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

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        false
    );

    assert!(result == true);
}

#[test]
fn test_variant_calling_8() {
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

    let result: bool = is_clusterable(
        &a,
        &b,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        false
    );

    assert!(result == true);
}

#[test]
fn test_variant_calling_9() {
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

    let result: bool = is_different(
        &a.graph_operation,
        &b.graph_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        false
    );

    assert!(result == false);
}

#[test]
fn test_variant_calling_10() {
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

    let result: bool = is_different(
        &a.graph_operation,
        &b.graph_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        false
    );

    assert!(result == false);
}

#[test]
fn test_variant_calling_11() {
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

    let result: bool = is_different(
        &a.graph_operation,
        &b.graph_operation,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        false
    );

    assert!(result == false);
}

#[test]
fn test_variant_calling_12() {
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

    let result: bool = is_different(
        &a.graph_operation,
        &b.graph_operation,
        0.05f32,
        0.95f32,
        2000,
        1000,
        1000,
        true,
        false
    );

    assert!(result == false);
}

#[test]
fn test_variant_calling_13() {
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

    let result: bool = is_different(
        &a.graph_operation,
        &b.graph_operation,
        0.05f32,
        0.95f32,
        2000,
        1000,
        1000,
        true,
        false
    );

    assert!(result == true);
}

#[test]
fn test_variant_calling_14() {
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
        1
    );

    assert!(variant_records_map.get(&(0,0,VariantType::Insertion, GraphOperationType::Downstream, GraphOperationType::Upstream)).unwrap().len() == 1);
    assert!(variant_records_map.get(&(1,1,VariantType::Insertion, GraphOperationType::Downstream, GraphOperationType::Upstream)).unwrap().len() == 1);
}

#[test]
fn test_variant_calling_15() {
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
    let variant_record_clusters: Vec<VariantRecordCluster> = sweep_clusters(
        &variant_records,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        2,
        false
    );

    assert!(variant_record_clusters.len() == 3);
}

#[test]
fn test_variant_calling_16() {
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
    let variant_record_clusters: Vec<VariantRecordCluster> = sweep_clusters(
        &variant_records,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        2,
        false
    );
    assert!(variant_record_clusters.len() == 1);
}

#[test]
fn test_variant_calling_17() {
    let mut a: Vec<VariantRecord> = Vec::new();
    let mut b: Vec<VariantRecord> = Vec::new();

    // DEL:chr1:1010001-1020000
    let go_a1: GraphOperation = GraphOperation::new(
        0,
        1010000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1020001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_a1: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_a1
    );
    a.push(vr_a1);

    // DEL:chr1:1005001-1020000
    let go_b1: GraphOperation = GraphOperation::new(
        0,
        1005000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1020001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_b1: VariantRecord = VariantRecord::new(
        2,
        2,
        3,
        go_b1
    );
    b.push(vr_b1);

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
        true,
        false
    );

    assert_eq!(a_only.is_empty(), true);
}
