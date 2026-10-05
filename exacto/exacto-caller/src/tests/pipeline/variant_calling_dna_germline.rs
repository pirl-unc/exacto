use csv::ReaderBuilder;
use std::fs;
use std::fs::File;
use std::path::Path;
use crate::io::builders::build_dna_variant_records;
use super::*;


#[test]
fn scga_mini_dna_001_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_002_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_003_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_004_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-004-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_005_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-005-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-005-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_006_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_007_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_008_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-008-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-008-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_009_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-009-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-009-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_010_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-010-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-010-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_011_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-011-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-011-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_012_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_013_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_014_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-014-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-014-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_015_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_016_identify_germline_dna_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam");
    let bam_bai_file = data_dir.join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam.bai");
    let fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file.to_str().unwrap(),
        bam_bai_file.to_str().unwrap(),
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
fn scga_mini_dna_001_identify_germline_dna_variants_with_region_of_two_chunks_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let mut options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    options.calling.chunk_size = 100_000;

    // The region is two chunks of 100,000 bases. The SNV at chr17:7,674,225 lies in the first chunk.
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file,
        bam_bai_file,
        fasta_file,
        &vec![("chr17", 7_600_001, 7_800_000)],
        &options,
        1,
        ""
    );
    let variant_call_set_no_regions: DNAVariantCallSet = identify_germline_dna_variants(
        bam_file,
        bam_bai_file,
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


/// Six reads whose first aligned base is chr18:1 and six whose last is chr18:10,000,000 (the last
/// base of the contig), each with a clip of 12 bases off the contig end. Such clips are
/// insertions at (0, 1) and (L, L + 1), and the run reads their depth at the base beside the
/// flank off the contig instead of stopping.
#[test]
fn identify_germline_dna_variants_calls_clips_on_the_first_and_last_base_of_a_contig() {
    use noodles_sam::alignment::io::Write;
    use tempfile::TempDir;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let clip: &str = "ACGTTGCAGTCA";
    let body: String = "ACGGTCATTGCA".repeat(17)[..200].to_string();
    let mut sam_text: String = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr17\tLN:10000000\n@SQ\tSN:chr18\tLN:10000000\n".to_string();
    for i in 0..6 {
        let flag: &str = if i % 2 == 0 { "0" } else { "16" };
        sam_text.push_str(&format!("first_{i}\t{flag}\tchr18\t1\t60\t12S200=\t*\t0\t0\t{clip}{body}\t{}\tcs:Z::200\n", "I".repeat(212)));
    }
    for i in 0..6 {
        let flag: &str = if i % 2 == 0 { "0" } else { "16" };
        sam_text.push_str(&format!("last_{i}\t{flag}\tchr18\t9999801\t60\t200=12S\t*\t0\t0\t{body}{clip}\t{}\tcs:Z::200\n", "I".repeat(212)));
    }
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("contig_edges.bam").to_str().unwrap().to_string();
    let bam_bai_file: String = format!("{bam_file}.bai");
    let mut sam_reader = noodles_sam::io::Reader::new(sam_text.as_bytes());
    let header: Header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&bam_bai_file, &bam::fs::index(&bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        &bam_file,
        &bam_bai_file,
        fasta_file,
        &vec![],
        &options,
        1,
        ""
    );

    let calls: Vec<(u32, u32, String, usize)> = variant_call_set
        .get_variant_calls()
        .iter()
        .map(|variant_call| {
            let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
            (operation.get_position_1(), operation.get_position_2(), operation.get_standardized_sequence(), variant_call.get_read_ids().len())
        })
        .collect();
    assert_eq!(
        calls,
        vec![(0, 1, "ACGTTGCAGTCA".to_string(), 6), (10_000_000, 10_000_001, "ACGTTGCAGTCA".to_string(), 6)]
    );
}


/// minimap2 writes supplementary records with hard clips unless given -Y. scga-mini-dna-007, whose
/// split reads hold the breakend chr17:3,491,600 -> chr17:6,085,001, rewritten that way gives the
/// same calls as with soft clips: the read sequence comes from the primary record, so a hard clip
/// stands for read bases as a soft clip does.
#[test]
fn scga_mini_dna_007_identify_germline_dna_variants_reads_hard_clipped_supplementary_records() {
    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let temp_dir = TempDir::new().unwrap();
    let hard_clipped_bam_file: String = temp_dir.path().join("hard_clipped.bam").to_str().unwrap().to_string();
    let hard_clipped_bam_bai_file: String = format!("{hard_clipped_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&hard_clipped_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    let mut num_hard_clipped: usize = 0;
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        if record.flags().is_supplementary() {
            let operations: Vec<Op> = record.cigar().as_ref().to_vec();
            let leading: usize = if operations[0].kind() == Kind::SoftClip { operations[0].len() } else { 0 };
            let trailing: usize = if operations[operations.len() - 1].kind() == Kind::SoftClip { operations[operations.len() - 1].len() } else { 0 };
            let length: usize = record.sequence().len();
            let sequence: Vec<u8> = record.sequence().as_ref()[leading..length - trailing].to_vec();
            let quality_scores: Vec<u8> = record.quality_scores().as_ref()[leading..length - trailing].to_vec();
            *record.sequence_mut() = sequence.into();
            *record.quality_scores_mut() = quality_scores.into();
            *record.cigar_mut() = operations
                .into_iter()
                .map(|op| if op.kind() == Kind::SoftClip { Op::new(Kind::HardClip, op.len()) } else { op })
                .collect();
            num_hard_clipped += 1;
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&hard_clipped_bam_bai_file, &bam::fs::index(&hard_clipped_bam_file).unwrap()).unwrap();
    assert_eq!(num_hard_clipped, 39);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let soft_clipped: DNAVariantCallSet = identify_germline_dna_variants(bam_file, bam_bai_file, fasta_file, &vec![], &options, 1, "");
    let hard_clipped: DNAVariantCallSet = identify_germline_dna_variants(&hard_clipped_bam_file, &hard_clipped_bam_bai_file, fasta_file, &vec![], &options, 1, "");

    let soft_clipped_rows: Vec<DNAVariantRecord> = build_dna_variant_records(&soft_clipped).collect();
    let hard_clipped_rows: Vec<DNAVariantRecord> = build_dna_variant_records(&hard_clipped).collect();
    assert_eq!(soft_clipped_rows.len(), 1);
    assert_eq!(hard_clipped_rows, soft_clipped_rows);
}


/// The index passed as its own path, under a name other than `<bam>.bai`: every pass reads it.
#[test]
fn scga_mini_dna_001_identify_germline_dna_variants_reads_an_index_stored_under_another_name() {
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let copied_bam_file: String = temp_dir.path().join("sample.bam").to_str().unwrap().to_string();
    let other_index_file: String = temp_dir.path().join("sample.index").to_str().unwrap().to_string();
    fs::copy(bam_file, &copied_bam_file).unwrap();
    fs::copy(bam_bai_file, &other_index_file).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let expected: DNAVariantCallSet = identify_germline_dna_variants(bam_file, bam_bai_file, fasta_file, &vec![], &options, 1, "");
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(&copied_bam_file, &other_index_file, fasta_file, &vec![], &options, 1, "");

    let expected_rows: Vec<DNAVariantRecord> = build_dna_variant_records(&expected).collect();
    let rows: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();
    assert_eq!(expected_rows.len(), 1);
    assert_eq!(rows, expected_rows);
}


/// A region subset of scga-mini-dna-007 (the records overlapping chr17:3,489,001-3,492,000). The
/// 20 supplementary records there belong to reads whose primary record lies past the breakend, at
/// chr17:6.07-6.09 Mb, so those reads keep only a supplementary record. The run goes through: a
/// read without its primary record has no read sequence and is left out.
#[test]
fn scga_mini_dna_007_identify_germline_dna_variants_runs_on_a_region_subset() {
    use noodles_sam::alignment::io::Write;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let subset_bam_file: String = temp_dir.path().join("subset.bam").to_str().unwrap().to_string();
    let subset_bam_bai_file: String = format!("{subset_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&subset_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    let mut num_records: usize = 0;
    let mut num_supplementary: usize = 0;
    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        let (Some(start), Some(end)) = (record.alignment_start(), record.alignment_end()) else {
            continue;
        };
        if record.reference_sequence_id() == Some(0) && start.get() <= 3_492_000 && end.get() >= 3_489_001 {
            num_records += 1;
            num_supplementary += usize::from(record.flags().is_supplementary());
            writer.write_alignment_record(&header, &record).unwrap();
        }
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&subset_bam_bai_file, &bam::fs::index(&subset_bam_file).unwrap()).unwrap();
    assert_eq!(num_records, 85);
    assert_eq!(num_supplementary, 20);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        &subset_bam_file, &subset_bam_bai_file, fasta_file, &vec![], &options, 1, ""
    );

    assert_eq!(variant_call_set.get_size(), 1);
}


/// A file with a header and no reads has nothing to call: an empty call set, and a table of no
/// rows.
#[test]
fn identify_germline_dna_variants_returns_no_call_for_a_file_without_reads() {
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let empty_bam_file: String = temp_dir.path().join("empty.bam").to_str().unwrap().to_string();
    let empty_bam_bai_file: String = format!("{empty_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&empty_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&empty_bam_bai_file, &bam::fs::index(&empty_bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
        &empty_bam_file, &empty_bam_bai_file, fasta_file, &vec![], &options, 1, ""
    );

    assert_eq!(variant_call_set.get_size(), 0);
    assert_eq!(build_dna_variant_records(&variant_call_set).count(), 0);
}


/// Mapping quality 255 means "not available". It is read as 255, as samtools reads it, so every
/// read passes the minimum and the calls are those of the file as written.
#[test]
fn scga_mini_dna_001_identify_germline_dna_variants_reads_mapping_quality_255() {
    use noodles_sam::alignment::io::Write;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let mapq_bam_file: String = temp_dir.path().join("mapq_255.bam").to_str().unwrap().to_string();
    let mapq_bam_bai_file: String = format!("{mapq_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&mapq_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        *record.mapping_quality_mut() = None;
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&mapq_bam_bai_file, &bam::fs::index(&mapq_bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    let expected: DNAVariantCallSet = identify_germline_dna_variants(bam_file, bam_bai_file, fasta_file, &vec![], &options, 1, "");
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(&mapq_bam_file, &mapq_bam_bai_file, fasta_file, &vec![], &options, 1, "");

    let expected_rows: Vec<DNAVariantRecord> = build_dna_variant_records(&expected).collect();
    let rows: Vec<DNAVariantRecord> = build_dna_variant_records(&variant_call_set).collect();
    assert_eq!(expected_rows.len(), 1);
    assert_eq!(rows, expected_rows);
}


/// A file aligned without --cs is refused before the first pass, with a message that names the
/// file and a read.
#[test]
#[should_panic(expected = "has no cs tag; align with minimap2 --cs.")]
fn scga_mini_dna_001_identify_germline_dna_variants_panics_for_a_file_without_cs_tags() {
    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::data::field::Tag;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let no_cs_bam_file: String = temp_dir.path().join("no_cs.bam").to_str().unwrap().to_string();
    let no_cs_bam_bai_file: String = format!("{no_cs_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&no_cs_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        record.data_mut().remove(&Tag::from([b'c', b's']));
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&no_cs_bam_bai_file, &bam::fs::index(&no_cs_bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    identify_germline_dna_variants(
        &no_cs_bam_file, &no_cs_bam_bai_file, fasta_file, &vec![], &options, 1, ""
    );
}


/// A file whose reads lie on a contig the reference genome lacks (here chr17 named 17, as
/// Ensembl names it) is refused before the first pass.
#[test]
#[should_panic(expected = "not in the reference genome FASTA file")]
fn scga_mini_dna_001_identify_germline_dna_variants_panics_for_a_contig_missing_from_the_reference() {
    use noodles_sam::alignment::io::Write;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let other_bam_file: String = temp_dir.path().join("other_contig.bam").to_str().unwrap().to_string();
    let other_bam_bai_file: String = format!("{other_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let other_header: Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:17\tLN:10000000\n@SQ\tSN:chr18\tLN:10000000\n"
        .parse()
        .unwrap();
    let mut writer = bam::io::Writer::new(File::create(&other_bam_file).unwrap());
    writer.write_header(&other_header).unwrap();
    for result in reader.record_bufs(&header) {
        writer.write_alignment_record(&other_header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&other_bam_bai_file, &bam::fs::index(&other_bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    identify_germline_dna_variants(
        &other_bam_file, &other_bam_bai_file, fasta_file, &vec![], &options, 1, ""
    );
}


/// A header that names a contig the reference genome lacks is refused even when the contig holds
/// no read (a genome-wide header over a reference of a few chromosomes): the check reads the
/// header, not the reads.
#[test]
#[should_panic(expected = "not in the reference genome FASTA file")]
fn scga_mini_dna_001_identify_germline_dna_variants_panics_for_a_contig_without_reads_missing_from_the_reference() {
    use noodles_sam::alignment::io::Write;
    use tempfile::TempDir;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let other_bam_file: String = temp_dir.path().join("empty_contig.bam").to_str().unwrap().to_string();
    let other_bam_bai_file: String = format!("{other_bam_file}.bai");
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let other_header: Header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr17\tLN:10000000\n@SQ\tSN:chr18\tLN:10000000\n@SQ\tSN:chrUn_KI270302v1\tLN:2274\n"
        .parse()
        .unwrap();
    let mut writer = bam::io::Writer::new(File::create(&other_bam_file).unwrap());
    writer.write_header(&other_header).unwrap();
    for result in reader.record_bufs(&header) {
        writer.write_alignment_record(&other_header, &result.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bai::fs::write(&other_bam_bai_file, &bam::fs::index(&other_bam_file).unwrap()).unwrap();

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
    identify_germline_dna_variants(
        &other_bam_file, &other_bam_bai_file, fasta_file, &vec![], &options, 1, ""
    );
}


/// Three runs of one input give the same table: the same ids, rows and read names in the same
/// order. scga-mini-dna-015 holds two calls, one of them a breakend held by split reads.
#[test]
fn scga_mini_dna_015_identify_germline_dna_variants_gives_the_same_table_on_three_runs() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let mut tables: Vec<Vec<DNAVariantRecord>> = Vec::new();
    for _ in 0..3 {
        let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
            bam_file, bam_bai_file, fasta_file, &vec![], &options, 1, ""
        );
        tables.push(build_dna_variant_records(&variant_call_set).collect());
    }

    assert_eq!(tables[0].len(), 2);
    assert_eq!(tables[1], tables[0]);
    assert_eq!(tables[2], tables[0]);
}


/// scga-mini-dna-015's breakend chr17:3,491,600 -> chr17:6,085,001 is held by 35 reads whose
/// records have position 1 from 3,491,589 to 3,492,224. A chunk size of 34,918 puts a chunk end
/// at 3,491,800, inside that spread. The window of that chunk owns the call and reads past the
/// chunk's end, so the call is made once and whole, as with the default chunk size.
#[test]
fn scga_mini_dna_015_identify_germline_dna_variants_does_not_depend_on_where_a_chunk_ends() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-015-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let mut variant_calls_by_chunk_size: Vec<Vec<String>> = Vec::new();
    for chunk_size in [100_000, 34_918] {
        let mut options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();
        options.calling.chunk_size = chunk_size;
        let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(
            bam_file,
            bam_bai_file,
            fasta_file,
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

    assert_eq!(variant_calls_by_chunk_size[0].len(), 2);
    assert!(variant_calls_by_chunk_size[0].iter().any(|variant_call| variant_call.starts_with("chr17:3491600:+:D:chr17:6085001:+:U:")));
    assert_eq!(variant_calls_by_chunk_size[1], variant_calls_by_chunk_size[0]);
}


/// A read flagged duplicate or QC-failed counts neither in the read depth nor in the read
/// support, so flagging reads gives the calls that removing them gives. The 31 reads of
/// scga-mini-dna-001 that carry its SNV chr17:7,674,225 C>A are flagged duplicate in one copy of
/// the file and left out of another. The two copies give the same calls, with the same depths,
/// and the SNV is not among them.
#[test]
fn scga_mini_dna_001_identify_germline_dna_variants_leaves_out_reads_flagged_duplicate() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::RecordBuf;
    use std::collections::HashSet;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    // The reads of the SNV.
    let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(bam_file, bam_bai_file, fasta_file, &vec![], &options, 1, "");
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
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
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
        let variant_call_set: DNAVariantCallSet = identify_germline_dna_variants(file, &bai_file, fasta_file, &vec![], &options, 1, "");
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
