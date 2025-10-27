use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_alignment_record_1() {
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let alignment_record_1: AlignmentRecord = AlignmentRecord::new(
        0,
        99,
        Strand::Forward,
        records_map.get(&1).unwrap().get(0).unwrap().clone()
    );

    let alignment_record_2: AlignmentRecord = alignment_record_1.clone();

    assert_eq!(alignment_record_1, alignment_record_2);
}
