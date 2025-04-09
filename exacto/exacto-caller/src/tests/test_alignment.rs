use bimap::BiMap;
use noodles_bam as bam;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_alignment_1() {
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

    let read_name: &str = "m64012_397004_551695/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_original_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_original_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        read_sequence,
        quality_scores,
        records_map.get(&read_id).unwrap().clone()
    );

    assert!(alignment.get_alignment_records_count() == 2);
}

#[test]
fn test_alignment_2() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_file: &str = reference_genome_fasta_path.to_str().unwrap();

    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );

    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );

    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);

    let mut alignments: Vec<Alignment> = Vec::new();
    for (read_id,records) in records_map.iter() {
        let original_read_sequence: Box<str> = get_original_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
        let original_base_quality_scores: Vec<u8> = get_original_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
        let alignment: Alignment = Alignment::new(
            *read_id,
            original_read_sequence,
            original_base_quality_scores,
            records.clone()
        );
        alignments.push(alignment);
    }

    let mut found: bool = false;
    for alignment in alignments.iter() {
        if alignment.read_id == *read_names_map.get_by_left(&"m64012_507476_774164/1/ccs".to_string().into_boxed_str()).unwrap() {
            let exons: Vec<TranscriptModelExon> = alignment.identify_exons(25);
            let splice_junctions: Vec<TranscriptModelSpliceJunction> = alignment.identify_splice_junctions(
                &chromosome_names_map,
                reference_genome_fasta_file,
                25
            );
            assert!(exons.len() == 11);
            assert!(splice_junctions.len() == 10);
            found = true;
        }
    }

    assert!(found);
}
