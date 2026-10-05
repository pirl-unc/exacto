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
use noodles_sam::alignment::io::Write as _;
use tempfile::TempDir;

use crate::prelude::*;


#[test]
fn test_bam_calculate_average_base_quality_score() {
    let base_quality_scores: Vec<u8> = vec![30,30,60,60];
    let average_base_quality_score: f32 = calculate_average_base_quality_score(&base_quality_scores);
    assert!(average_base_quality_score == 45.0);
}

#[test]
fn test_bam_check_bam_end_of_file() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    check_bam_end_of_file(bam_file);

    // The same file cut at 30% of its bytes, as an interrupted copy leaves it. Before the check,
    // every reader returned the records before the cut without a message.
    let temp_dir = TempDir::new().unwrap();
    let cut_bam_path = temp_dir.path().join("cut.bam");
    let cut_bam_file: &str = cut_bam_path.to_str().unwrap();
    let bytes: Vec<u8> = fs::read(bam_file).unwrap();
    fs::write(cut_bam_file, &bytes[..bytes.len() * 3 / 10]).unwrap();

    let read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    let results = [
        std::panic::catch_unwind(|| { index_bam_records(cut_bam_file, true, 1); }),
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            fetch_all_bam_records(cut_bam_file, bam_bai_file, &read_names_map, 1);
        })),
        std::panic::catch_unwind(|| { get_read_names_passing_mapping_quality(cut_bam_file, bam_bai_file, 1, 20); }),
        std::panic::catch_unwind(|| { get_read_names(cut_bam_file, bam_bai_file, 1); }),
        std::panic::catch_unwind(|| { get_bam_depths_map(cut_bam_file, &format!("{}.bai", cut_bam_file), 1); }),
        std::panic::catch_unwind(|| { get_chromosome_lengths(cut_bam_file); }),
    ];
    for result in results {
        let error = result.unwrap_err();
        let message: &String = error.downcast_ref::<String>().unwrap();
        assert!(message.contains("does not end with the BGZF end-of-file block"), "{message}");
    }
}

#[test]
fn test_bam_create_chromosome_names_map_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosomes_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    assert!(chromosomes_map.len() == 2);
    assert!(chromosomes_map.contains_left("chr17") == true);
    assert!(chromosomes_map.contains_left("chr18") == true);
}

#[test]
fn test_bam_create_chromosome_names_map_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    assert_eq!(chromosome_names_map.len(), 2);
}

#[test]
fn test_bam_create_read_names_map() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    assert!(read_names_map.len() == 68);
}

#[test]
fn test_bam_fetch_all_bam_records() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(bam_file, bam_bai_file, &read_names_map, 1);
    assert!(read_names_map.len() == 68);
    assert!(records.keys().len() == 68);
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
        assert!(records.len() == 1);
        assert!(read_names.contains(&read_name) == true);
    }
}

#[test]
fn test_bam_fetch_bam_records_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);

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

    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        7_600_000,
        7_700_000,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );
    assert!(read_names_map.len() == 68);
    assert!(records.keys().len() == 68);
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(read_id).unwrap().clone();
        assert!(records.len() == 1);
        assert!(read_names.contains(&read_name) == true);
    }
}

#[test]
fn test_bam_fetch_bam_records_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

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

    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        7_600_000,
        7_700_000,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );
    // dna-006 holds an inversion of chr17:7,701,201-7,717,000; 32 of the 57 reads in the region
    // are split by it into a primary and a supplementary record. Its one unmapped read is left out.
    assert!(read_names_map.len() == 164);
    assert!(records.keys().len() == 57);
    assert_eq!(records.values().filter(|records| records.len() == 2).count(), 32);
    assert_eq!(records.values().filter(|records| records.len() == 1).count(), 25);
}

#[test]
fn test_bam_fetch_bam_records_3() {
    // A BAM cut to a region keeps the secondary and supplementary records whose primary record
    // lies outside the region. r_secondary_only and r_supplementary_only have only those;
    // r_primary has a secondary record at chrS:1801.
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("region_subset.bam").to_str().unwrap().to_string();
    let bam_bai_file: String = format!("{bam_file}.bai");
    let sam: String = [
        "@HD\tVN:1.6\tSO:coordinate".to_string(),
        "@SQ\tSN:chrS\tLN:2000".to_string(),
        format!("r_primary\t0\tchrS\t101\t60\t100M\t*\t0\t0\t{}\t*", "A".repeat(100)),
        format!("r_chimeric\t0\tchrS\t201\t60\t50M50S\t*\t0\t0\t{}{}\t*", "C".repeat(50), "G".repeat(50)),
        format!("r_secondary_only\t256\tchrS\t601\t0\t100M\t*\t0\t0\t*\t*"),
        format!("r_supplementary_only\t2048\tchrS\t1001\t60\t50S50M\t*\t0\t0\t{}\t*", "T".repeat(100)),
        format!("r_chimeric\t2048\tchrS\t1201\t60\t50H50M\t*\t0\t0\t{}\t*", "G".repeat(50)),
        format!("r_primary\t256\tchrS\t1801\t0\t100M\t*\t0\t0\t*\t*"),
    ].join("\n") + "\n";
    let mut sam_reader = noodles_sam::io::Reader::new(sam.as_bytes());
    let header: Header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&bam_bai_file, &bam::fs::index(&bam_file).unwrap()).unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(&bam_file, true, 1);
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(&bam_bai_file).unwrap();

    // The region of the secondary-only and supplementary-only reads: nothing to fetch, no panic.
    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader, &header, &index, "chrS".into(), 550, 1100, &record_positions_map, &read_names_map, 7, 1
    );
    assert!(records.is_empty());

    // A secondary record brings no read into its region.
    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader, &header, &index, "chrS".into(), 1750, 1950, &record_positions_map, &read_names_map, 7, 1
    );
    assert!(records.is_empty());

    // The whole contig: every read with a primary record, each with its primary and supplementary
    // records, and every one gives its read as sequenced.
    let records: HashMap<usize, Vec<Record>> = fetch_bam_records(
        &mut reader, &header, &index, "chrS".into(), 1, 2000, &record_positions_map, &read_names_map, 7, 1
    );
    let mut records_per_read: Vec<(Box<str>, usize)> = records
        .iter()
        .map(|(read_id, records)| (read_names_map.get_by_right(read_id).unwrap().clone(), records.len()))
        .collect();
    records_per_read.sort();
    assert_eq!(records_per_read, vec![("r_chimeric".into(), 2), ("r_primary".into(), 1)]);
    for records in records.values() {
        assert_eq!(get_bam_fastx_read_sequence(records).len(), 100);
    }
}

#[test]
fn test_bam_generate_buffered_regions() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_buffered_regions(
        bam_file,
        &vec!["chr17".into(), "chr18".into()],
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
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_regions(
        &vec!["chr17".into(), "chr18".into()],
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
#[should_panic(expected = "chunk_size must be > 0")]
fn test_bam_generate_regions_chunk_size_0() {
    // A chunk size of 0 stops at once; it used to repeat one region for ever.
    let chromosome_lengths: HashMap<Box<str>, u32> = HashMap::from([("chr1".into(), 1_000)]);
    generate_regions(&vec!["chr1".into()], &chromosome_lengths, 0);
}

#[test]
#[should_panic(expected = "chunk_size must be > 0")]
fn test_bam_generate_buffered_regions_chunk_size_0() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    generate_buffered_regions(bam_file, &vec!["chr17".into()], 0, 10_000);
}

#[test]
fn test_bam_get_alignment_end_position() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let alignment_end: u32 = get_alignment_end_position(record);
    assert!(alignment_end == 7686329);
}

#[test]
fn test_bam_get_alignment_mapping_quality() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let mapping_quality: u16 = get_alignment_mapping_quality(record);
    assert!(mapping_quality == 60);
}

#[test]
fn test_bam_get_alignment_mapping_quality_255() {
    // 255 means "not available" (STAR writes it for every unique mapper); noodles reads it as None.
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("mapq.bam").to_str().unwrap().to_string();
    let bam_bai_file: String = format!("{bam_file}.bai");
    let sam: String = [
        "@HD\tVN:1.6\tSO:coordinate".to_string(),
        "@SQ\tSN:chrS\tLN:2000".to_string(),
        format!("r_mapq60\t0\tchrS\t101\t60\t100M\t*\t0\t0\t{}\t*", "A".repeat(100)),
        format!("r_mapq10\t0\tchrS\t201\t10\t100M\t*\t0\t0\t{}\t*", "C".repeat(100)),
        format!("r_mapq255\t0\tchrS\t301\t255\t100M\t*\t0\t0\t{}\t*", "G".repeat(100)),
    ].join("\n") + "\n";
    let mut sam_reader = noodles_sam::io::Reader::new(sam.as_bytes());
    let header: Header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&bam_bai_file, &bam::fs::index(&bam_file).unwrap()).unwrap();

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam_file)
        .unwrap();
    reader.read_header().unwrap();
    let mapping_qualities: Vec<u16> = reader
        .records()
        .map(|result| get_alignment_mapping_quality(&result.unwrap()))
        .collect();
    assert_eq!(mapping_qualities, vec![60, 10, 255]);

    let read_names: HashSet<Box<str>> = get_read_names_passing_mapping_quality(&bam_file, &bam_bai_file, 1, 20);
    let expected: HashSet<Box<str>> = ["r_mapq60".into(), "r_mapq255".into()].into_iter().collect();
    assert_eq!(read_names, expected);
}

#[test]
fn test_bam_get_aligned_sequence_from_cigar() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(record);
    // The read's 15,879 bases less its 4 soft-clipped ones.
    assert!(aligned_sequence.len() == 15875);
}

#[test]
fn test_bam_get_aligned_sequence_from_cigar_hard_clip() {
    // minimap2 hard-clips supplementary records unless run with -Y.
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("hard_clip.bam").to_str().unwrap().to_string();
    let sam: String = [
        "@HD\tVN:1.6\tSO:coordinate".to_string(),
        "@SQ\tSN:chrS\tLN:2000".to_string(),
        format!("r_chimeric\t0\tchrS\t201\t60\t50M50S\t*\t0\t0\t{}{}\t*", "C".repeat(50), "G".repeat(50)),
        format!("r_chimeric\t2048\tchrS\t1201\t60\t50H50M\t*\t0\t0\t{}\t*", "G".repeat(50)),
        format!("r_chimeric\t2064\tchrS\t1501\t60\t20H30M50H\t*\t0\t0\t{}{}\t*", "C".repeat(20), "A".repeat(10)),
    ].join("\n") + "\n";
    let mut sam_reader = noodles_sam::io::Reader::new(sam.as_bytes());
    let header: Header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam_file)
        .unwrap();
    reader.read_header().unwrap();
    let aligned_sequences: Vec<Box<str>> = reader
        .records()
        .map(|result| get_aligned_sequence_from_cigar(&result.unwrap()))
        .collect();
    let expected: Vec<Box<str>> = vec![
        "C".repeat(50).into(),
        "G".repeat(50).into(),
        format!("{}{}", "T".repeat(10), "G".repeat(20)).into(),
    ];
    assert_eq!(aligned_sequences, expected);
}

#[test]
fn test_bam_get_alignment_start_position() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let alignment_start: u32 = get_alignment_start_position(record);
    assert!(alignment_start == 7670402);
}

#[test]
fn test_bam_get_alignment_strand() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let strand: Strand = get_alignment_strand(record);
    assert!(strand == Strand::Forward);
}

#[test]
fn test_bam_get_bam_header() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let header: Header = get_bam_header(bam_file);
    assert!(header.reference_sequences().len() == 2);
}

#[test]
fn test_bam_get_chromosome_lengths() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    assert!(chromosome_lengths.keys().len() == 2);
}

#[test]
fn test_bam_get_chromosome_names() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    assert!(chromosome_names.len() == 2);
}

#[test]
fn test_bam_get_cigar_operations() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &bam::Record = records.get(&read_id).unwrap().get(0).unwrap();
    let cigar_ops: Vec<(Kind, u32)> = get_cigar_operations(record);
    assert!(cigar_ops.len() == 177);
    assert!(cigar_ops[0].0 == Kind::SoftClip);
    assert!(cigar_ops[0].1 == 4);
    assert!(cigar_ops[1].0 == Kind::SequenceMatch);
    assert!(cigar_ops[1].1 == 144);
    assert!(cigar_ops[2].0 == Kind::Deletion);
    assert!(cigar_ops[2].1 == 1);
}

#[test]
fn test_bam_get_fastx_base_quality_scores() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: Record = records.get(&read_id).unwrap().get(0).unwrap().clone();
    let scores: Vec<u8> = get_bam_fastx_base_quality_scores(&vec![record]);
    assert!(scores.len() == 15879);
}

#[test]
fn test_bam_get_fastx_read_sequence() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    assert!(sequence.len() == 15879);
}

#[test]
fn test_bam_get_left_softclipping_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-007-tumor_scga-mini-dna-007-tumor-3_1/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(get_left_softclipping(record).0 == false);
}

#[test]
fn test_bam_get_left_softclipping_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-007-tumor_scga-mini-dna-007-tumor-3_1/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(1).unwrap();
    assert!(get_left_softclipping(record).0 == true);
}

#[test]
fn test_bam_get_right_softclipping_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-007-tumor_scga-mini-dna-007-tumor-3_1/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(get_right_softclipping(record).0 == true);
}

#[test]
fn test_bam_get_right_softclipping_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-007-tumor_scga-mini-dna-007-tumor-3_1/2/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(1).unwrap();
    assert!(get_right_softclipping(record).0 == false);
}

#[test]
fn test_bam_get_primary_alignment_base_quality_scores() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let record_refs: Vec<&Record> = records.iter().collect();
    let scores: Vec<u8> = get_primary_alignment_base_quality_scores(&record_refs);
    assert!(scores.len() == 15879);
}

#[test]
fn test_bam_get_primary_alignment_read_sequence() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let records: &Vec<Record> = records.get(&read_id).unwrap();
    let record_refs: Vec<&Record> = records.iter().collect();
    let sequence: Box<str> = get_primary_alignment_read_sequence(&record_refs);
    assert!(sequence.len() == 15879);
}

#[test]
fn test_bam_get_read_names() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names(bam_file, bam_bai_file, 1);
    assert!(read_names.len() == 68);
}

#[test]
fn test_bam_get_read_names_passing_mapping_quality_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names: HashSet<Box<str>> = get_read_names_passing_mapping_quality(
        bam_file,
        bam_bai_file,
        1,
        60
    );
    assert!(read_names.len() == 68);
}

#[test]
fn test_bam_get_read_names_passing_mapping_quality_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
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
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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
    assert!(read_names.len() == 358);
}

#[test]
fn test_bam_get_read_sequence() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let sequence: Box<str> = get_read_sequence(record);
    assert!(sequence.len() == 15879);
}

#[test]
fn test_bam_get_tag_value() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-001-tumor_chunk_0000/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    let cs_tag: Option<Box<str>> = get_tag_value(record, "cs");
    assert!(cs_tag.is_some());
    assert!(&*cs_tag.unwrap() == ":190+c:226-a:268-g:302+c:261+c:21~ct918ac:107~ct2819ac:74~ct92ac:137~ct343ac:44*ca:65~ct568ac:113~ct81ac:184~ct757ac:279~ct109ac:22~ct117ac:102~ct10754ac:51");
}

#[test]
fn test_bam_has_soft_clipping_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-001-tumor_chunk_0000/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_soft_clipping(record) == false);
}

#[test]
fn test_bam_has_soft_clipping_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-004-tumor_scga-mini-dna-004-tumor-2_1/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_soft_clipping(record) == true);
}

#[test]
fn test_bam_has_splicing_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-004-tumor_scga-mini-dna-004-tumor-2_1/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_splicing(record) == false);
}

#[test]
fn test_bam_has_splicing_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-001-tumor_chunk_0000/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_splicing(record) == true);
}

#[test]
fn test_bam_has_tag() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-001-tumor_chunk_0000/1/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(has_tag(record, "cs") == true);
}

#[test]
fn test_bam_index_bam_records_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );
    assert!(record_positions_map.len() == 68);
    assert!(read_names_map.len() == 68);
}

#[test]
fn test_bam_index_bam_records_2() {
    // Reads with no primary record in the file are left out; the others keep all their primary
    // and supplementary records, and no secondary one.
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("region_subset.bam").to_str().unwrap().to_string();
    let sam: String = [
        "@HD\tVN:1.6\tSO:coordinate".to_string(),
        "@SQ\tSN:chrS\tLN:2000".to_string(),
        format!("r_primary\t0\tchrS\t101\t60\t100M\t*\t0\t0\t{}\t*", "A".repeat(100)),
        format!("r_chimeric\t0\tchrS\t201\t60\t50M50S\t*\t0\t0\t{}{}\t*", "C".repeat(50), "G".repeat(50)),
        format!("r_secondary_only\t256\tchrS\t601\t0\t100M\t*\t0\t0\t*\t*"),
        format!("r_supplementary_only\t2048\tchrS\t1001\t60\t50S50M\t*\t0\t0\t{}\t*", "T".repeat(100)),
        format!("r_chimeric\t2048\tchrS\t1201\t60\t50H50M\t*\t0\t0\t{}\t*", "G".repeat(50)),
        format!("r_primary\t256\tchrS\t1801\t0\t100M\t*\t0\t0\t*\t*"),
    ].join("\n") + "\n";
    let mut sam_reader = noodles_sam::io::Reader::new(sam.as_bytes());
    let header: Header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(fs::File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);

    let (record_positions_map, read_names_map) = index_bam_records(&bam_file, true, 1);
    let mut records_per_read: Vec<(Box<str>, usize)> = read_names_map
        .iter()
        .map(|(read_name, read_id)| (read_name.clone(), record_positions_map.get(read_id).unwrap().len()))
        .collect();
    records_per_read.sort();
    assert_eq!(records_per_read, vec![("r_chimeric".into(), 2), ("r_primary".into(), 1)]);
    assert_eq!(record_positions_map.len(), 2);
}


#[test]
fn test_bam_is_aligned_to_reverse_strand_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/25/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(is_aligned_to_reverse_strand(record) == false);
}

#[test]
fn test_bam_is_aligned_to_reverse_strand_2() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-2_1/14/ccs").unwrap();
    let record: &Record = records.get(&read_id).unwrap().get(0).unwrap();
    assert!(is_aligned_to_reverse_strand(record) == true);
}

#[test]
fn test_bam_get_bam_depths_map() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let depths_map: HashMap<Box<str>, Vec<u32>> = get_bam_depths_map(bam_file, &format!("{}.bai", bam_file), 2);
    let depth: u32 = *depths_map.get("chr17").unwrap().get(7_673_000 - 1).unwrap();
    assert!(depth == 66);
}

#[test]
fn test_bam_read_depths_get_strands() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let positions_map: HashMap<Box<str>, Vec<u32>> = HashMap::from([("chr17".into(), vec![7_673_000])]);
    let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, bam_bai_file, &positions_map, 1_000);
    let (fwd_count, rev_count): (u32, u32) = read_depths.get_strands("chr17", 7_673_000);
    assert!(fwd_count == 36);
    assert!(rev_count == 30);
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
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );
    let read_id: usize = *read_names_map.get_by_left("scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-2_1/14/ccs").unwrap();
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

#[test]
fn test_bam_read_depths_match_the_whole_contig_maps() {
    // dna-004 holds a 120 bp insertion (long I operations), dna-006 an inversion (supplementary records).
    for name in ["dna-004-tumor", "dna-006-tumor"] {
        let bam_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(format!("alignment/scga-mini-{name}_minimap2_sorted.bam"))).unwrap();
        let bam_file: &str = bam_path.to_str().unwrap();
        let bam_bai_file: String = format!("{bam_file}.bai");
        let depths_map: HashMap<Box<str>, Vec<u32>> = get_bam_depths_map(bam_file, &bam_bai_file, 2);

        // Every covered position with 100 bases either side, and the two positions off each contig.
        let mut positions: Vec<(Box<str>, u32)> = Vec::new();
        for (contig, depths) in depths_map.iter() {
            let covered: Vec<usize> = (0..depths.len())
                .filter(|&i| depths[i] > 0)
                .collect();
            if let (Some(&first), Some(&last)) = (covered.first(), covered.last()) {
                let start: u32 = first.saturating_sub(100) as u32 + 1;
                let end: u32 = (last + 100).min(depths.len() - 1) as u32 + 1;
                positions.extend((start..=end).map(|position| (contig.clone(), position)));
            }
            positions.push((contig.clone(), 0));
            positions.push((contig.clone(), depths.len() as u32 + 1));
        }
        assert!(positions.len() > 10_000);

        let mut positions_map: HashMap<Box<str>, Vec<u32>> = HashMap::new();
        for (contig, position) in positions.iter() {
            positions_map.entry(contig.clone()).or_default().push(*position);
        }
        let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &bam_bai_file, &positions_map, 1_000);

        for (contig, position) in positions {
            let length: usize = depths_map[&contig].len();
            let offset: usize = (position as usize).clamp(1, length) - 1;
            assert_eq!(read_depths.get_depth(&contig, position), depths_map[&contig][offset], "{name} {contig}:{position}");
        }
        assert_eq!(
            get_bam_max_depth(bam_file, &bam_bai_file, 2),
            depths_map.values().flat_map(|depths| depths.iter().copied()).max().unwrap()
        );
    }
}

#[test]
fn test_bam_read_depths_strands_add_up_to_depths_map() {
    // Depth and strand counts are taken over the same reads (supplementary kept; secondary,
    // QC-failed and duplicate left out; deleted bases not counted), so at every base the forward
    // and the reverse count add up to the depth.
    for name in ["dna-001-tumor", "dna-004-tumor", "dna-006-tumor"] {
        let bam_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(format!("alignment/scga-mini-{name}_minimap2_sorted.bam"))).unwrap();
        let bam_file: &str = bam_path.to_str().unwrap();
        let bam_bai_file: String = format!("{bam_file}.bai");
        let depths_map: HashMap<Box<str>, Vec<u32>> = get_bam_depths_map(bam_file, &bam_bai_file, 2);

        // Every base from 100 before the first covered base to 100 after the last.
        let mut positions_map: HashMap<Box<str>, Vec<u32>> = HashMap::new();
        for (contig, depths) in depths_map.iter() {
            let covered: Vec<usize> = (0..depths.len())
                .filter(|&i| depths[i] > 0)
                .collect();
            if let (Some(&first), Some(&last)) = (covered.first(), covered.last()) {
                let start: u32 = first.saturating_sub(100) as u32 + 1;
                let end: u32 = (last + 100).min(depths.len() - 1) as u32 + 1;
                positions_map.insert(contig.clone(), (start..=end).collect());
            }
        }
        let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &bam_bai_file, &positions_map, 1_000);

        let mut num_covered: usize = 0;
        for (contig, positions) in positions_map.iter() {
            for &position in positions {
                let (forward, reverse): (u32, u32) = read_depths.get_strands(contig, position);
                let depth: u32 = depths_map[contig][position as usize - 1];
                assert_eq!(forward + reverse, depth, "{name} {contig}:{position}");
                num_covered += (depth > 0) as usize;
            }
        }
        assert!(num_covered > 1_000, "{name}");
    }
}
