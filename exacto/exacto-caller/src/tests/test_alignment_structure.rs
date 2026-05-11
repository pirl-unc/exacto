use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use crate::prelude::*;


#[test]
fn test_alignment_structure_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut found: bool = false;
    for base in alignment.get_alignment_structure().get_bases() {
        if *base.get_kind() == AlignmentStructureBaseKind::Mismatch {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_382982_262550/2/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut found: bool = false;
    for base in alignment.get_alignment_structure().get_bases() {
        if *base.get_kind() == AlignmentStructureBaseKind::Insertion {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_3() {
    let bam_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_478275_464661/2/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut found: bool = false;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Deletion {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_4() {
    let bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_767230_904257/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut num_breakpoint: usize = 0;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Breakpoint {
            num_breakpoint += 1;
        }
    }

    assert_eq!(num_breakpoint, 2);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_5() {
    let bam_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_283345_480209/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut num_breakpoint: usize = 0;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Breakpoint {
            num_breakpoint += 1;
        }
    }

    for ((r1,r2), event) in alignment_structure.get_events().iter() {
        if *event.get_kind() == AlignmentStructureEventKind::Breakpoint {
            let b1 = alignment_structure.get_base(*r1);
            let b2 = alignment_structure.get_base(*r2);
            assert!(b1.get_reference_position().unwrap() == 4637155);
            assert!(b2.get_reference_position().unwrap() == 7674880);
        }
    }

    assert_eq!(num_breakpoint, 1);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_6() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_825713_352116/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut num_breakpoints: usize = 0;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Breakpoint {
            num_breakpoints += 1;
        }
    }

    assert_eq!(num_breakpoints, 1);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_7() {
    let bam_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_291012_248279/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut num_breakpoints: usize = 0;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Breakpoint {
            num_breakpoints += 1;
        }
    }

    assert_eq!(num_breakpoints, 1);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), false);
}

#[test]
fn test_alignment_structure_8() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_507476_774164/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut num_mismatch: usize = 0;
    for base in alignment.get_alignment_structure().get_bases() {
        if *base.get_kind() == AlignmentStructureBaseKind::Mismatch {
            num_mismatch += 1;
        }
    }

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    assert_eq!(num_mismatch, 1);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_9() {
    let bam_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_822724_603243/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut num_insertion: usize = 0;
    for base in alignment.get_alignment_structure().get_bases() {
        if *base.get_kind() == AlignmentStructureBaseKind::Insertion {
            num_insertion += 1;
        }
    }

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    assert_eq!(num_insertion, 10);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_10() {
    let bam_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
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

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_264855_304921/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let mut num_deletion: usize = 0;
    for event in alignment.get_alignment_structure().get_events().values() {
        if event.get_kind() == &AlignmentStructureEventKind::Deletion {
            num_deletion += 1;
        }
    }

    assert_eq!(num_deletion, 1);
    assert_eq!(alignment.get_alignment_structure().is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_11() {
    let bam_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

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

    let read_name: &str = "m64012_535544_475898/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000698746.1").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_mismatch: usize = 0;
    for base in alignment_structure.get_bases() {
        if *base.get_kind() == AlignmentStructureBaseKind::Mismatch {
            num_mismatch += 1;
        }
    }

    let mut num_fusion_gene: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::FusionGene {
            num_fusion_gene += 1;
        }
    }

    let variant_records = alignment_structure.identify_variant_records(
        30,
        30,
        AnalyteType::RNA
    );

    let mut fusion_gene_exists: bool = false;
    for variant_record in variant_records.iter() {
        if *variant_record.get_variant_type() == VariantType::FusionGene {
            fusion_gene_exists = true;
            break;
        }
    }

    assert_eq!(num_mismatch, 1);
    assert_eq!(num_fusion_gene, 1);
    assert_eq!(fusion_gene_exists, true);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 18);
    assert_eq!(introns.len(), 17);
}

#[test]
fn test_alignment_structure_12() {
    let bam_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

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

    let read_name: &str = "m64012_561742_839878/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_ref_bases_skipped: usize = 0;
    for event in alignment_structure.get_events().values() {
        for chromsome_id in event.get_skipped_reference_bases().keys() {
            num_ref_bases_skipped += event
                .get_skipped_reference_bases()
                .get(chromsome_id)
                .unwrap()
                .len();
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::NonCanonicalSplicing {
            num_noncanonical_splicing += 1;
        }
    }

    assert_eq!(num_ref_bases_skipped, 106);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_13() {
    let bam_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_124525_407996/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_ref_bases_skipped: usize = 0;
    for event in alignment_structure.get_events().values() {
        for chromsome_id in event.get_skipped_reference_bases().keys() {
            num_ref_bases_skipped += event
                .get_skipped_reference_bases()
                .get(chromsome_id)
                .unwrap()
                .len();
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::NonCanonicalSplicing {
            num_noncanonical_splicing += 1;
        }
    }

    assert_eq!(num_ref_bases_skipped, 72);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_14() {
    let bam_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

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

    let read_name: &str = "m64012_924107_174289/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_ref_bases_skipped: usize = 0;
    for event in alignment_structure.get_events().values() {
        for chromsome_id in event.get_skipped_reference_bases().keys() {
            num_ref_bases_skipped += event
                .get_skipped_reference_bases()
                .get(chromsome_id)
                .unwrap()
                .len();
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::NonCanonicalSplicing {
            num_noncanonical_splicing += 1;
        }
    }

    assert_eq!(num_ref_bases_skipped, 279);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 10);
    assert_eq!(introns.len(), 9);
}

#[test]
fn test_alignment_structure_15() {
    let bam_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_924958_759981/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_ref_bases_cryptic: usize = 0;
    for base in alignment_structure.get_bases() {
        if *base.get_context().as_ref().unwrap() == AlignmentStructureBaseContext::Intronic {
            num_ref_bases_cryptic += 1;
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::NonCanonicalSplicing {
            num_noncanonical_splicing += 1;
        }
    }

    assert_eq!(num_ref_bases_cryptic, 41);
    assert_eq!(num_noncanonical_splicing, 2);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 12);
    assert_eq!(introns.len(), 11);
}

#[test]
fn test_alignment_structure_16() {
    let bam_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

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

    let read_name: &str = "m64012_721712_133913/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_ref_bases_intron: usize = 0;
    for base in alignment_structure.get_bases() {
        if *base.get_context().as_ref().unwrap() == AlignmentStructureBaseContext::Intronic {
            num_ref_bases_intron += 1;
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::NonCanonicalSplicing {
            num_noncanonical_splicing += 1;
        }
    }

    assert_eq!(num_ref_bases_intron, 18);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(introns.len(), 10);
}

#[test]
fn test_alignment_structure_17() {
    let bam_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_288476_571946/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_1: &Transcript = gene_annotator.get_transcript("ENST00000263087.9").unwrap();
    let reference_transcript_2: &Transcript = gene_annotator.get_transcript("ENST00000570791.5").unwrap();
    let reference_transcript_3: &Transcript = gene_annotator.get_transcript("ENST00000333813.4").unwrap();
    let reference_transcript_sequence_1: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_1,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    let reference_transcript_sequence_2: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_2,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    let reference_transcript_sequence_3: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_3,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![
            &reference_transcript_sequence_1,
            &reference_transcript_sequence_2,
            &reference_transcript_sequence_3
        ],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_fusion_gene: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::FusionGene {
            num_fusion_gene += 1;
        }
    }

    let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(
        30,
        30,
        AnalyteType::RNA
    );

    assert_eq!(num_fusion_gene, 2);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 19);
    assert_eq!(introns.len(), 16);

    assert_eq!(variant_records.len(), 35);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 1295600);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 3801187);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(0).unwrap().get_sequence(), "G"); // overlapping alignment
    assert_eq!(variant_records.get(29).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(29).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(29).unwrap().get_position_1(), 3761101);
    assert_eq!(variant_records.get(29).unwrap().get_position_2(), 7727201);
    assert_eq!(variant_records.get(29).unwrap().get_operation_1(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(29).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(29).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(29).unwrap().get_sequence(), "G"); // overlapping alignment
}

#[test]
fn test_alignment_structure_18() {
    let bam_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();


    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_175366_924183/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let exons: Vec<TranscriptModelExon> = alignment_structure.identify_exons("");
    let introns: Vec<TranscriptModelIntron> = alignment_structure.identify_introns(
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_1: &Transcript = gene_annotator.get_transcript("ENST00000263092.11").unwrap();
    let reference_transcript_2: &Transcript = gene_annotator.get_transcript("ENST00000250113.12").unwrap();
    let reference_transcript_3: &Transcript = gene_annotator.get_transcript("ENST00000355530.7").unwrap();
    let reference_transcript_sequence_1: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_1,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    let reference_transcript_sequence_2: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_2,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    let reference_transcript_sequence_3: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript_3,
        &chromosome_names_map,
        reference_genome_fasta_file
    );
    alignment_structure.contextualize(
        "",
        &vec![
            &reference_transcript_sequence_1,
            &reference_transcript_sequence_2,
            &reference_transcript_sequence_3,
        ],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_fusion_gene: usize = 0;
    for event in alignment_structure.get_events().values() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::FusionGene {
            num_fusion_gene += 1;
        }
    }

    let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(
        30,
        30,
        AnalyteType::RNA
    );

    assert_eq!(num_fusion_gene, 2);
    assert_eq!(alignment_structure.is_spliced(), true);
    assert_eq!(exons.len(), 27);
    assert_eq!(introns.len(), 24);

    assert_eq!(variant_records.len(), 17);
    assert_eq!(variant_records.get(4).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_position_1(), 2464208);
    assert_eq!(variant_records.get(4).unwrap().get_position_2(), 4433940);
    assert_eq!(variant_records.get(4).unwrap().get_operation_1(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(4).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(4).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(4).unwrap().get_sequence(), "");
    assert_eq!(variant_records.get(6).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_position_1(), 4453100);
    assert_eq!(variant_records.get(6).unwrap().get_position_2(), 7603799);
    assert_eq!(variant_records.get(6).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(6).unwrap().get_operation_2(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(6).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(6).unwrap().get_sequence(), "");
}

#[test]
fn test_alignment_structure_19() {
    let bam_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
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

    let read_name: &str = "m64012_324970_273886/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect()
    );

    let mut alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000254719.10").unwrap();
    let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
        reference_transcript,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    alignment_structure.contextualize(
        "",
        &vec![&reference_transcript_sequence],
        &gene_annotator,
        &chromosome_names_map
    );

    let mut num_backsplicing: usize = 0;
    for ((i,j),event) in alignment_structure.get_events() {
        if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::BackSplicing {
            num_backsplicing += 1;
        }
    }

    assert_eq!(num_backsplicing, 1);
    assert_eq!(alignment_structure.is_spliced(), true);
}
