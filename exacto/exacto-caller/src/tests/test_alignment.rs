use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_alignment_1() {
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
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert_eq!(alignment.get_alignment_records_count(), 1);
    assert_eq!(alignment.get_base_quality_scores().len(), 19070);
    assert_eq!(alignment.get_read_id(), read_id);
    assert_eq!(alignment.get_read_length(), alignment.get_base_quality_scores().len());
    assert_eq!(alignment.get_read_sequence().len(), alignment.get_base_quality_scores().len());
}

#[test]
fn test_alignment_2() {
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

    let read_name: &str = "m64012_382982_262550/2/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
}

#[test]
fn test_alignment_3() {
    let bam_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_478275_464661/2/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
}

#[test]
fn test_alignment_4() {
    let bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_767230_904257/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert_eq!(alignment.get_alignment_records_count(), 3);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_alignment_records().get(2).unwrap().reference_strand, Strand::Forward);
}

#[test]
fn test_alignment_5() {
    let bam_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_283345_480209/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 2);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Forward);
}

#[test]
fn test_alignment_6() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_825713_352116/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 2);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Forward);
}

#[test]
fn test_alignment_7() {
    let bam_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_291012_248279/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 2);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_8() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_507476_774164/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_9() {
    let bam_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_822724_603243/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_10() {
    let bam_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_264855_304921/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_11() {
    let bam_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_535544_475898/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
}

#[test]
fn test_alignment_12() {
    let bam_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_561742_839878/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_13() {
    let bam_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_124525_407996/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_14() {
    let bam_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_924107_174289/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_15() {
    let bam_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_924958_759981/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_16() {
    let bam_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_721712_133913/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 1);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_17() {
    let bam_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_288476_571946/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 3);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_alignment_records().get(2).unwrap().reference_strand, Strand::Forward);
}

#[test]
fn test_alignment_18() {
    let bam_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_175366_924183/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 3);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(2).unwrap().reference_strand, Strand::Reverse);
}

#[test]
fn test_alignment_19() {
    let bam_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam.bai");
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

    let read_name: &str = "m64012_324970_273886/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    assert!(alignment.get_alignment_records_count() == 2);
    assert_eq!(alignment.get_alignment_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_alignment_records().get(1).unwrap().reference_strand, Strand::Forward);
}
