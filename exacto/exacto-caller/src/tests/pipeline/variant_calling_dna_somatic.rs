use csv::ReaderBuilder;
use std::fs;
use std::fs::File;
use std::path::Path;
use crate::io::builders::build_dna_variant_records;
use super::*;


#[test]
fn scga_mini_dna_001_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-001-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_002_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-002-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-002-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-002-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_003_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-003-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-003-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-003-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_004_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-004-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-004-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-004-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_005_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-005-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-005-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-005-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-005-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-005-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_006_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-006-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-006-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-006-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_007_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-007-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-007-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-007-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_008_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-008-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-008-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-008-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-008-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-008-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_009_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-009-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-009-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-009-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-009-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-009-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_010_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-010-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-010-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-010-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-010-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-010-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_011_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-011-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-011-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-011-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-011-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-011-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_012_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-012-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-012-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes and operations. The two copies meet across a G that both sides hold (7679901 and 7673300), so the junction
    // can be written one base either way; the caller writes the G as inserted sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-012-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && (r.position_1 as i64 - row[2].parse::<i64>().unwrap()).abs() <= 1
                && (r.position_2 as i64 - row[6].parse::<i64>().unwrap()).abs() <= 1
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_013_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-013-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-013-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-013-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_014_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-014-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-014-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-014-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-014-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-014-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_015_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-015-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-015-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-015-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_016_identify_somatic_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-016-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-016-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    let records: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();

    // Compare against the ground truth: each row is exactly one call, with the same
    // chromosomes, operations, positions and sequence.
    let tsv_file = data_dir.join("simulation/ground_truth/scga-mini-dna-016-tumor_ground_truth.tsv");
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(tsv_file).unwrap();
    let mut num_rows: usize = 0;
    for result in reader.records() {
        let row = result.unwrap();
        num_rows += 1;
        let num_matches: usize = records
            .iter()
            .filter(|r| {
                &*r.chromosome_1 == &row[1]
                && &*r.operation_1 == &row[4]
                && &*r.chromosome_2 == &row[5]
                && &*r.operation_2 == &row[8]
                && r.position_1 == row[2].parse::<u32>().unwrap()
                && r.position_2 == row[6].parse::<u32>().unwrap()
                && r.sequence.to_uppercase() == row[11].to_uppercase()
            })
            .count();
        assert_eq!(num_matches, 1, "ground truth row {:?} must match exactly one call", row);
    }
    assert_eq!(records.len(), num_rows);
}

#[test]
fn scga_mini_dna_001_identify_somatic_dna_variants_with_region_of_two_chunks_returns_matches() {
    let tumor_bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam");
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam.bai");
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let mut options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    options.calling.chunk_size = 100_000;

    // The region is two chunks of 100,000 bases. The SNV at chr17:7,674,225 lies in the first chunk.
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        vec![normal_bam_file],
        vec![normal_bam_bai_file],
        fasta_file,
        &vec![("chr17", 7_600_001, 7_800_000)],
        &options,
        1,
        ""
    );
    let variant_call_set_no_regions: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        vec![normal_bam_file],
        vec![normal_bam_bai_file],
        fasta_file,
        &vec![],
        &options,
        1,
        ""
    );

    let mut variant_calls: Vec<String> = Vec::new();
    for variant_call in variant_call_set.get_variant_calls().iter() {
        let mut read_ids = variant_call.get_read_ids();
        read_ids.sort();
        variant_calls.push(format!(
            "{} {:?}",
            variant_call.get_consensus_graph_operation().as_named_boxed_str(&variant_call_set.chromosome_names_map),
            read_ids
        ));
    }
    variant_calls.sort();
    let mut variant_calls_no_regions: Vec<String> = Vec::new();
    for variant_call in variant_call_set_no_regions.get_variant_calls().iter() {
        let mut read_ids = variant_call.get_read_ids();
        read_ids.sort();
        variant_calls_no_regions.push(format!(
            "{} {:?}",
            variant_call.get_consensus_graph_operation().as_named_boxed_str(&variant_call_set_no_regions.chromosome_names_map),
            read_ids
        ));
    }
    variant_calls_no_regions.sort();

    assert!(variant_calls.len() == 1);
    assert!(variant_calls[0].starts_with("chr17:7674224:+:D:chr17:7674226:+:U:A:1:SNV "));
    assert_eq!(variant_calls, variant_calls_no_regions);
}

/// A tumor that is its own control has no somatic variant.
///
/// scga-mini-dna-007 holds a breakend joining chr17:3,491,600 to chr17:6,085,001, which its reads
/// hold as split alignments. Against the normal the tumor calls it.
#[test]
fn scga_mini_dna_007_identify_somatic_dna_variants_returns_no_call_when_control_is_the_tumor() {
    let tumor_bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let mut options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    options.filtering.min_reads = 3;
    options.filtering.min_total_depth = 3;
    options.error_model.sequencing_error = 0.001;
    options.error_model.slippage_prob = 0.002;

    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        vec![tumor_bam_file],
        vec![tumor_bam_bai_file],
        fasta_file,
        &vec![],
        &options,
        1,
        ""
    );

    assert!(variant_call_set.get_size() == 0);
}

/// A tumor that is its own control has no somatic variant.
///
/// scga-mini-dna-004 holds an insertion of 120 bases after chr17:7,674,224. Each of its 32 carrier
/// reads spans it and holds it as one insertion, of 117 to 121 bases; the control's insertions
/// are the same alleles. Against the normal the tumor calls it.
#[test]
fn scga_mini_dna_004_identify_somatic_dna_variants_returns_no_call_when_control_is_the_tumor() {
    let tumor_bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam.bai");
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let mut options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    options.filtering.min_reads = 3;
    options.filtering.min_total_depth = 3;
    options.error_model.sequencing_error = 0.001;
    options.error_model.slippage_prob = 0.002;

    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        vec![tumor_bam_file],
        vec![tumor_bam_bai_file],
        fasta_file,
        &vec![],
        &options,
        1,
        ""
    );

    assert!(variant_call_set.get_size() == 0);
}

/// A control that holds the tumor's variant as insertions only still holds it.
///
/// scga-mini-dna-013 holds an inverted duplication of 6,600 bases. A read that spans the
/// duplicated copy holds it as one insertion; a read that ends inside it holds it as split
/// alignments. The tumor has both kinds, so its insertions are retyped as breakpoints and the
/// calls are breakpoints. The control here is made of the tumor's reads of the first kind.
#[test]
fn scga_mini_dna_013_identify_somatic_dna_variants_returns_no_call_when_control_holds_insertions_only() {
    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use tempfile::TempDir;

    let tumor_bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    // Write the control: the tumor's reads with one record and an insertion of 100 bases or more.
    let temp_dir = TempDir::new().unwrap();
    let control_bam_file: String = temp_dir.path().join("control.bam").to_str().unwrap().to_string();
    let control_bam_bai_file: String = format!("{control_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(tumor_bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let records: Vec<bam::Record> = reader.records().map(|result| result.unwrap()).collect();
    let mut num_records: HashMap<String, usize> = HashMap::new();
    for record in records.iter() {
        *num_records.entry(record.name().unwrap().to_string()).or_insert(0) += 1;
    }
    let mut writer = bam::io::Writer::new(File::create(&control_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    let mut num_control_reads: usize = 0;
    for record in records.iter() {
        let has_long_insertion: bool = record
            .cigar()
            .iter()
            .map(|result| result.unwrap())
            .any(|op| op.kind() == Kind::Insertion && op.len() >= 100);
        if has_long_insertion && num_records[&record.name().unwrap().to_string()] == 1 {
            writer.write_alignment_record(&header, record).unwrap();
            num_control_reads += 1;
        }
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&control_bam_bai_file, &bam::fs::index(&control_bam_file).unwrap()).unwrap();
    assert_eq!(num_control_reads, 20);

    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        vec![control_bam_file.as_str()],
        vec![control_bam_bai_file.as_str()],
        fasta_file,
        &vec![("chr17", 7_600_001, 7_700_000)],
        &options,
        1,
        ""
    );

    assert!(variant_call_set.get_size() == 0);
}


/// `max_control_reads` is the most control reads a variant may have and still be called
/// somatic. The control holds two of the tumour's reads that carry its SNV, both read from one
/// strand: with 2 the call stays, with 1 it is subtracted. Only the SNV is counted: the control
/// lacks the normal's reads, so a 1-base deletion at 7,679,230 that the normal subtracts is
/// called against it either way.
#[test]
fn scga_mini_dna_001_identify_somatic_dna_variants_allows_max_control_reads_in_the_control() {
    use noodles_sam::alignment::io::Write;
    use std::collections::HashSet;
    use tempfile::TempDir;

    let tumor_bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam")).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: String = format!("{tumor_bam_file}.bai");
    let normal_bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam")).unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: String = format!("{normal_bam_file}.bai");
    let fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    // The somatic call against the normal, and two of its reads.
    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file, &tumor_bam_bai_file, vec![normal_bam_file], vec![normal_bam_bai_file.as_str()],
        fasta_file, &vec![], &options, 1, ""
    );
    assert_eq!(variant_call_set.get_size(), 1);
    let mut forward_read_ids: Vec<usize> = variant_call_set.get_variant_calls()[0]
        .get_variant_records()
        .iter()
        .filter(|variant_record| *variant_record.get_strand_1() == Strand::Forward)
        .map(|variant_record| variant_record.get_read_id())
        .collect();
    forward_read_ids.sort();
    let control_read_names: HashSet<String> = forward_read_ids
        .iter()
        .take(2)
        .map(|read_id| variant_call_set.read_names_map.get_by_right(read_id).unwrap().to_string())
        .collect();
    assert_eq!(control_read_names.len(), 2);

    // The control: every record of those two reads.
    let temp_dir = TempDir::new().unwrap();
    let control_bam_file: String = temp_dir.path().join("control.bam").to_str().unwrap().to_string();
    let control_bam_bai_file: String = format!("{control_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(tumor_bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&control_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.records() {
        let record: bam::Record = result.unwrap();
        if control_read_names.contains(&record.name().unwrap().to_string()) {
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&control_bam_bai_file, &bam::fs::index(&control_bam_file).unwrap()).unwrap();

    let mut num_calls: Vec<usize> = Vec::new();
    for max_control_reads in [2, 1] {
        let mut options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
        options.subtraction.max_control_reads = max_control_reads;
        let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
            tumor_bam_file, &tumor_bam_bai_file, vec![control_bam_file.as_str()], vec![control_bam_bai_file.as_str()],
            fasta_file, &vec![], &options, 1, ""
        );
        num_calls.push(
            variant_call_set
                .get_variant_calls()
                .iter()
                .filter(|variant_call| variant_call.get_consensus_graph_operation().get_position_1() == 7_674_224)
                .count()
        );
    }

    assert_eq!(num_calls, vec![1, 0]);
}


/// scga-mini-dna-015's breakend chr17:3,491,600 -> chr17:6,085,001 is held by reads whose records
/// have position 1 from 3,491,589 to 3,492,224. A chunk size of 34,918 puts a chunk end at
/// 3,491,800, inside that spread. The window of that chunk owns the call and reads past the
/// chunk's end, so the somatic call is made once and whole, as with the default chunk size.
#[test]
fn scga_mini_dna_015_identify_somatic_dna_variants_does_not_depend_on_where_a_chunk_ends() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-015-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-015-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let mut variant_calls_by_chunk_size: Vec<Vec<String>> = Vec::new();
    for chunk_size in [100_000, 34_918] {
        let mut options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();
        options.calling.chunk_size = chunk_size;
        let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
            tumor_bam_file.to_str().unwrap(),
            tumor_bam_bai_file.to_str().unwrap(),
            vec![normal_bam_file.to_str().unwrap()],
            vec![normal_bam_bai_file.to_str().unwrap()],
            fasta_file.to_str().unwrap(),
            &vec![],
            &options,
            1,
            ""
        );
        let mut variant_calls: Vec<String> = variant_call_set
            .get_variant_calls()
            .iter()
            .map(|variant_call| format!(
                "{} {:?}",
                variant_call.get_consensus_graph_operation().as_named_boxed_str(&variant_call_set.chromosome_names_map),
                variant_call.get_read_ids()
            ))
            .collect();
        variant_calls.sort();
        variant_calls_by_chunk_size.push(variant_calls);
    }

    assert!(variant_calls_by_chunk_size[0].iter().any(|variant_call| variant_call.starts_with("chr17:3491600:+:D:chr17:6085001:+:U:")));
    assert_eq!(variant_calls_by_chunk_size[1], variant_calls_by_chunk_size[0]);
}


/// A tumor read flagged duplicate counts neither in the read depth nor in the read support, so
/// flagging reads gives the calls that removing them gives. The 31 reads of scga-mini-dna-001
/// that carry its SNV chr17:7,674,225 C>A are flagged duplicate in one copy of the tumor file and
/// left out of another. Against the normal, the two copies give the same calls, and the SNV is
/// not among them.
#[test]
fn scga_mini_dna_001_identify_somatic_dna_variants_leaves_out_reads_flagged_duplicate() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::RecordBuf;
    use std::collections::HashSet;

    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let tumor_bam_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let tumor_bam_bai_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let normal_bam_file = data_dir.join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam");
    let normal_bam_bai_file = data_dir.join("alignment/scga-mini-dna-001-normal_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let options: IdentifySomaticDNAVariantsOptions = IdentifySomaticDNAVariantsOptions::default();

    // The reads of the SNV.
    let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
        tumor_bam_file.to_str().unwrap(),
        tumor_bam_bai_file.to_str().unwrap(),
        vec![normal_bam_file.to_str().unwrap()],
        vec![normal_bam_bai_file.to_str().unwrap()],
        fasta_file.to_str().unwrap(),
        &vec![],
        &options,
        1,
        ""
    );
    assert_eq!(variant_call_set.get_size(), 1);
    let alternate_read_names: HashSet<String> = variant_call_set.get_variant_calls()[0]
        .get_read_ids()
        .iter()
        .map(|read_id| variant_call_set.read_names_map.get_by_right(read_id).unwrap().to_string())
        .collect();
    assert_eq!(alternate_read_names.len(), 31);

    // One copy with every record of those reads flagged duplicate, one without the reads.
    let directory: tempfile::TempDir = tempfile::tempdir().unwrap();
    let flagged_bam_file: String = directory.path().join("flagged.bam").to_str().unwrap().to_string();
    let removed_bam_file: String = directory.path().join("removed.bam").to_str().unwrap().to_string();
    let mut reader = bam::io::reader::Builder::default().build_from_path(&tumor_bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut flagged_writer = bam::io::Writer::new(File::create(&flagged_bam_file).unwrap());
    let mut removed_writer = bam::io::Writer::new(File::create(&removed_bam_file).unwrap());
    flagged_writer.write_header(&header).unwrap();
    removed_writer.write_header(&header).unwrap();
    for result in reader.records() {
        let record: bam::Record = result.unwrap();
        let mut record_buf: RecordBuf = RecordBuf::try_from_alignment_record(&header, &record).unwrap();
        if alternate_read_names.contains(&record.name().unwrap().to_string()) {
            *record_buf.flags_mut() = record_buf.flags() | Flags::DUPLICATE;
        } else {
            removed_writer.write_alignment_record(&header, &record_buf).unwrap();
        }
        flagged_writer.write_alignment_record(&header, &record_buf).unwrap();
    }
    flagged_writer.try_finish().unwrap();
    removed_writer.try_finish().unwrap();
    drop(flagged_writer);
    drop(removed_writer);

    let mut variant_calls_by_file: Vec<Vec<String>> = Vec::new();
    for file in [&flagged_bam_file, &removed_bam_file] {
        let bai_file: String = format!("{file}.bai");
        bai::fs::write(&bai_file, &bam::fs::index(file).unwrap()).unwrap();
        let variant_call_set: DNAVariantCallSet = identify_somatic_dna_variants(
            file,
            &bai_file,
            vec![normal_bam_file.to_str().unwrap()],
            vec![normal_bam_bai_file.to_str().unwrap()],
            fasta_file.to_str().unwrap(),
            &vec![],
            &options,
            1,
            ""
        );
        let mut variant_calls: Vec<String> = variant_call_set
            .get_variant_calls()
            .iter()
            .map(|variant_call| {
                let read_names: Vec<String> = variant_call
                    .get_read_ids()
                    .iter()
                    .map(|read_id| variant_call_set.read_names_map.get_by_right(read_id).unwrap().to_string())
                    .collect();
                assert!(read_names.iter().all(|read_name| !alternate_read_names.contains(read_name)));
                format!(
                    "{} depth {} {:?}",
                    variant_call.get_consensus_graph_operation().as_named_boxed_str(&variant_call_set.chromosome_names_map),
                    variant_call.get_total_depth(),
                    read_names
                )
            })
            .collect();
        variant_calls.sort();
        variant_calls_by_file.push(variant_calls);
    }

    assert_eq!(variant_calls_by_file[0], variant_calls_by_file[1]);
    assert!(variant_calls_by_file[0].iter().all(|variant_call| !variant_call.starts_with("chr17:7674224:")));
}
