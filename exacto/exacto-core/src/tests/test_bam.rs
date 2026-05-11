use bimap::BiMap;
use noodles_bam as bam;
use noodles_bam::{bai, Record};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::Header;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use noodles_bam::bai::Index;
use tempfile::TempDir;

use crate::prelude::*;


#[test]
fn test_bam_calculate_average_base_quality_score() {
    let base_quality_scores: Vec<u8> = vec![30,30,60,60];
    let average_base_quality_score: f32 = calculate_average_base_quality_score(&base_quality_scores);
    assert!(average_base_quality_score == 45.0);
}

#[test]
fn test_bam_create_chromosome_names_map_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosomes_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    assert!(chromosomes_map.len() == 2);
    assert!(chromosomes_map.contains_left("chr17") == true);
    assert!(chromosomes_map.contains_left("chr18") == true);
}

#[test]
fn test_bam_create_chromosome_names_map_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    assert_eq!(chromosome_names_map.len(), 2);
}

#[test]
fn test_bam_create_read_names_map() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    assert!(read_names_map.len() == 6);
}

#[test]
fn test_bam_fetch_all_bam_records() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(bam_file, bam_bai_file, &read_names_map, 1);
    assert!(read_names_map.len() == 6);
    assert!(records.keys().len() == 6);
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
        assert!(records.len() == 1);
        assert!(read_names.contains(&read_name) == true);
    }
}

#[test]
fn test_bam_fetch_bam_records_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17",
        7_600_000,
        7_700_000,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );
    assert!(read_names_map.len() == 6);
    assert!(records.keys().len() == 6);
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
        assert!(records.len() == 1);
        assert!(read_names.contains(&read_name) == true);
    }
}

#[test]
fn test_bam_fetch_bam_records_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17",
        7_600_000,
        7_700_000,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );
    assert!(read_names_map.len() == 9);
    assert!(records.keys().len() == 6);
    let read_names_softclipped: HashSet<Box<str>> = [
        "m64012_825713_352116/1/ccs",
        "m64012_400530_589417/2/ccs",
        "m64012_925457_739219/3/ccs",
    ]
    .into_iter()
    .map(Into::into)
    .collect();
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
        if read_names_softclipped.contains(&read_name) == true {
            assert!(records.len() == 2);
        } else {
            assert!(records.len() == 1);
        }
    }
}

#[test]
fn test_bam_generate_buffered_regions() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_buffered_regions(
        bam_file,
        &vec!["chr17", "chr18"],
        1_000_000,
        10_000
    );
    assert!(regions.len() == 2);
    assert!(regions.contains_key("chr17") == true);
    assert!(regions.contains_key("chr18") == true);
    assert!(regions.get("chr17").unwrap().len() == 10);
    assert!(regions.get("chr18").unwrap().len() == 10);
}

#[test]
fn test_bam_generate_regions() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_regions(
        &vec!["chr17", "chr18"],
        &chromosome_lengths,
        1_000_000,
    );
    assert!(regions.len() == 2);
    assert!(regions.contains_key("chr17") == true);
    assert!(regions.contains_key("chr18") == true);
    assert!(regions.get("chr17").unwrap().len() == 10);
    assert!(regions.get("chr18").unwrap().len() == 10);
}

#[test]
fn test_bam_get_alignment_end_position() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let alignment_end: u32 = get_alignment_end_position(record);
    assert!(alignment_end == 7687490);
}

#[test]
fn test_bam_get_alignment_mapping_quality() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let mapping_quality: u16 = get_alignment_mapping_quality(record);
    assert!(mapping_quality == 60);
}

#[test]
fn test_bam_get_aligned_sequence_from_cigar() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(record);
    assert!(aligned_sequence.len() == 19070);
}

#[test]
fn test_bam_get_alignment_start_position() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let alignment_start: u32 = get_alignment_start_position(record);
    assert!(alignment_start == 7668421);
}

#[test]
fn test_bam_get_alignment_strand() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let strand: Strand = get_alignment_strand(record);
    assert!(strand == Strand::Forward);
}

#[test]
fn test_bam_get_bam_header() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let header: Header = get_bam_header(bam_file);
    assert!(header.reference_sequences().len() == 2);
}

#[test]
fn test_bam_get_chromosome_lengths() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    assert!(chromosome_lengths.keys().len() == 2);
}

#[test]
fn test_bam_get_chromosome_names() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    assert!(chromosome_names.len() == 2);
}

#[test]
fn test_bam_get_cigar_operations() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let cigar_ops: Vec<(Kind, u32)> = get_cigar_operations(record);
    assert!(cigar_ops.len() == 3);
    assert!(cigar_ops[0].0 == Kind::SequenceMatch);
    assert!(cigar_ops[0].1 == 5804);
    assert!(cigar_ops[1].0 == Kind::SequenceMismatch);
    assert!(cigar_ops[1].1 == 1);
    assert!(cigar_ops[2].0 == Kind::SequenceMatch);
    assert!(cigar_ops[2].1 == 13265);
}

#[test]
fn test_bam_get_fastx_base_quality_scores() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: Record = records.get(&read_id).unwrap().get(0).unwrap().clone();
    let scores: Vec<u8> = get_fastx_base_quality_scores(&vec![record]);
    assert!(scores.len() == 19070);
}

#[test]
fn test_bam_get_fastx_read_sequence() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let sequence: Box<str> = get_fastx_read_sequence(&records);
    assert!(sequence.len() == 19070);
}

#[test]
fn test_bam_get_left_softclipping_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_825713_352116/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(get_left_softclipping(record).0 == false);
}

#[test]
fn test_bam_get_left_softclipping_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_825713_352116/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(1).unwrap();
    assert!(get_left_softclipping(record).0 == true);
}

#[test]
fn test_bam_get_right_softclipping_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_825713_352116/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(get_right_softclipping(record).0 == true);
}

#[test]
fn test_bam_get_right_softclipping_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_825713_352116/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(1).unwrap();
    assert!(get_right_softclipping(record).0 == false);
}

#[test]
fn test_bam_get_primary_alignment_base_quality_scores() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let record_refs: Vec<&Record> = records.iter().collect();
    let scores: Vec<u8> = get_primary_alignment_base_quality_scores(&record_refs);
    assert!(scores.len() == 19070);
}

#[test]
fn test_bam_get_primary_alignment_read_sequence() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let record_refs: Vec<&Record> = records.iter().collect();
    let sequence: Box<str> = get_primary_alignment_read_sequence(&record_refs);
    assert!(sequence.len() == 19070);
}

#[test]
fn test_bam_get_read_names() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);
    assert!(read_names.len() == 6);
}

#[test]
fn test_bam_get_read_names_passing_mapping_quality_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names_passing_mapping_quality(
        bam_file,
        bam_bai_file,
        1,
        60
    );
    assert!(read_names.len() == 6);
}

#[test]
fn test_bam_get_read_names_passing_mapping_quality_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names_passing_mapping_quality(
        bam_file,
        bam_bai_file,
        1,
        100
    );
    assert!(read_names.len() == 0);
}

#[test]
fn test_bam_get_read_names_with_splicing_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let gencode: Gencode = Gencode::new(
        gtf_file,
        "hg38",
        "v41",
        Some(vec!["protein_coding"].into_iter().collect()),
        Some(vec![1,2].into_iter().collect()),
        Some(vec!["protein_coding"].into_iter().collect()),
        Some(vec![1,2].into_iter().collect())
    );
    let read_names: HashSet<Box<str>> = get_read_names_with_splicing(
        bam_file,
        bam_bai_file,
        &gencode,
        1
    );
    assert!(read_names.len() == 0);
}

#[test]
fn test_bam_get_read_names_with_splicing_2() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let gencode: Gencode = Gencode::new(
        gtf_file,
        "hg38",
        "v41",
        Some(vec!["protein_coding"].into_iter().collect()),
        Some(vec![1,2].into_iter().collect()),
        Some(vec!["protein_coding"].into_iter().collect()),
        Some(vec![1,2].into_iter().collect())
    );
    let read_names: HashSet<Box<str>> = get_read_names_with_splicing(
        bam_file,
        bam_bai_file,
        &gencode,
        1
    );
    assert!(read_names.len() == 2);
}

#[test]
fn test_bam_get_read_sequence() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let sequence: Box<str> = get_read_sequence(record);
    assert!(sequence.len() == 19070);
}

#[test]
fn test_bam_get_tag_value() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_507476_774164/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let cs_tag: Option<Box<str>> = get_tag_value(record, "cs");
    assert!(cs_tag.is_some());
    assert!(&*cs_tag.unwrap() == ":1270~ct918ac:107~ct2819ac:74~ct92ac:137~ct343ac:44*ca:65~ct568ac:113~ct81ac:184~ct757ac:279~ct109ac:22~ct117ac:102~ct10754ac:114");
}

#[test]
fn test_bam_has_soft_clipping_1() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_507476_774164/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_soft_clipping(record) == false);
}

#[test]
fn test_bam_has_soft_clipping_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_767230_904257/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_soft_clipping(record) == true);
}

#[test]
fn test_bam_has_splicing_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_767230_904257/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_splicing(record) == false);
}

#[test]
fn test_bam_has_splicing_2() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_507476_774164/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_splicing(record) == true);
}

#[test]
fn test_bam_has_tag() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_507476_774164/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_tag(record, "cs") == true);
}

#[test]
fn test_bam_index_bam_records_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );
    assert!(record_positions_map.len() == 6);
    assert!(read_names_map.len() == 6);
}


#[test]
fn test_bam_is_aligned_to_reverse_strand_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_325382_158010/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(is_aligned_to_reverse_strand(record) == false);
}

#[test]
fn test_bam_is_aligned_to_reverse_strand_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_899988_730625/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(is_aligned_to_reverse_strand(record) == true);
}

#[test]
fn test_bam_get_bam_depths_map() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let depths_map: HashMap<Box<str>, Vec<u32>> = get_bam_depths_map(bam_file, 2);
    let depth: u32 = *depths_map.get("chr17").unwrap().get(7_673_000 - 1).unwrap();
    assert!(depth == 6);
}

#[test]
fn test_bam_get_bam_strands_map() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let strands_map: HashMap<Box<str>, Vec<(u32, u32)>> = get_bam_strands_map(bam_file, 2);
    let fwd_count: u32 = strands_map.get("chr17").unwrap().get(7_673_000 - 1).unwrap().0;
    let rev_count: u32 = strands_map.get("chr17").unwrap().get(7_673_000 - 1).unwrap().1;
    assert!(fwd_count == 3);
    assert!(rev_count == 3);
}

#[test]
fn test_bam_split_regions_1() {
    let regions: Vec<(&str, u32, u32)> = vec![("chr17", 1_000_001, 2_000_000)];
    let regions_split = split_regions(&regions, 100_000);
    assert!(regions_split.len() == 10);
}

#[test]
fn test_bam_split_regions_2() {
    let regions: Vec<(&str, u32, u32)> = vec![("chr17", 1_000_000, 2_000_000)];
    let regions_split = split_regions(&regions, 100_000);
    assert!(regions_split.len() == 11);
}

#[test]
fn test_bam_split_regions_3() {
    let regions: Vec<(&str, u32, u32)> = vec![("chr17", 100_000, 200_000)];
    let regions_split = split_regions(&regions, 1_000_000);
    assert!(regions_split.len() == 1);
}

#[test]
fn test_bam_write_bam_file() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("m64012_899988_730625/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();

    let temp_dir = TempDir::new().unwrap();
    let output_bam_file: String = temp_dir.path().join("test.bam").to_str().unwrap().to_string();
    let output_bai_file: String = temp_dir.path().join("test.bam.bai").to_str().unwrap().to_string();

    let header: Header = get_bam_header(bam_file);
    
    write_bam_file(
        output_bam_file.as_str(),
        output_bai_file.as_str(),
        &header,
        &vec![record]
    );
}
