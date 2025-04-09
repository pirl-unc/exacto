use bimap::BiMap;
use bam::io::Reader;
use noodles_bam as bam;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::common::bam::*;


#[test]
fn test_fetch_all_bam_records_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_ids_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_ids_map,
        1
    );
    assert!(records_map.keys().len() == 6);
}

#[test]
fn test_fetch_bam_records_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_ids_map: BiMap<Box<str>,usize> = create_read_names_map(
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
        &read_ids_map,
        1
    );
    assert!(records_map.keys().len() == 6);
}

#[test]
fn test_get_aligned_sequence_1() {
    let path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let full_path = fs::canonicalize(path).unwrap();
    let file = File::open(full_path.to_str().unwrap()).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        let primary_alignment_read_sequence: Box<str> = get_primary_alignment_read_sequence(&vec![&record]);
        let aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(&record);
        assert_eq!(aligned_sequence.len(), 19070);
    } else {
        panic!("Failed to read the first BAM record.");
    }
}

#[test]
fn test_get_alignment_end_position_1() {
    let path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let full_path = fs::canonicalize(path).unwrap();
    let file = File::open(full_path.to_str().unwrap()).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        let alignment_end: usize = get_alignment_end_position(&record);
        assert_eq!(alignment_end, 7687490);
    } else {
        panic!("Failed to read the first BAM record.");
    }
}

#[test]
fn test_get_alignment_start_position_1() {
    let path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let full_path = fs::canonicalize(path).unwrap();
    let file = File::open(full_path.to_str().unwrap()).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        let alignment_end: usize = get_alignment_start_position(&record);
        assert_eq!(alignment_end, 7668421);
    } else {
        panic!("Failed to read the first BAM record.");
    }
}

#[test]
fn test_get_chromosomes_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);
    let chr17: u16 = *chromosome_names_map.get_by_left("chr17").unwrap();
    let chr18: u16 = *chromosome_names_map.get_by_left("chr18").unwrap();
    assert!(chr17 == 0 || chr17 == 1);
    assert!(chr18 == 0 || chr18 == 1);
}

#[test]
fn test_get_chromosome_ids_names_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);
    assert_eq!(chromosome_names_map.len(), 2);
}

#[test]
fn test_get_chromosome_lengths_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosomes: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    assert_eq!(chromosomes.keys().len(), 2);
}

#[test]
fn test_get_original_read_sequence() {
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
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_original_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    assert!(read_sequence != "".into());
}

#[test]
fn test_get_read_sequence_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let file = File::open(bam_file).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        let read_sequence: Box<str> = get_read_sequence(&record);
        assert_eq!(read_sequence.len(), 19070);
    } else {
        panic!("Failed to read the first BAM record.");
    }
}

#[test]
fn test_get_tag_value_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let file = File::open(bam_file).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        if let Some(value) = get_tag_value(&record, "cs") {
            let cs_tag: String = value.to_string();
            assert!(cs_tag == ":5804*ca:13265");
        } else {
            panic!("CS tag not found.");
        }
    } else {
        panic!("Failed to read the first BAM record.");
    }
}

#[test]
fn test_is_aligned_to_reverse_strand_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let file = File::open(bam_file).unwrap();
    let mut reader = Reader::new(BufReader::new(file));
    let _ = reader.read_header();

    // Read the first record
    if let Some(Ok(record)) = reader.records().next() {
        assert!(is_aligned_to_reverse_strand(&record) == false);
    } else {
        panic!("Failed to read the first BAM record.");
    }
}
