use csv::ReaderBuilder;
use polars::prelude::*;
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::prelude::*;


#[test]
fn tset_variant_calling_dna_1() {
    let bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        3,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-001-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("variant_sequence").str().to_uppercase().eq(lit(variant_sequence.to_uppercase())))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_2() {
    let bam_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        3,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-002-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("variant_sequence").str().to_uppercase().eq(lit(variant_sequence.to_uppercase())))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_3() {
    let bam_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        3,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-003-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_4() {
    let bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        3,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 2);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-004-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 1000)))
            .filter(col("position_1").lt_eq(lit(position_1 + 1000)))
            .filter(col("position_2").gt_eq(lit(position_2 - 1000)))
            .filter(col("position_2").lt_eq(lit(position_2 + 1000)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_5() {
    let bam_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        3,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-005-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_6() {
    let bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        1,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17","chr18"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-006-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_7() {
    let bam_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam");
    let bam_bai_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let variant_call_set: DNAVariantCallSet = identify_dna_variants(
        bam_file,
        bam_bai_file,
        1,
        25,
        25,
        0.5f32,
        0.5f32,
        20000,
        1000,
        1000,
        1,
        vec!["chr17","chr18"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-007-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_8() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-001-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-001-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_9() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-002-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-002-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-002-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-002-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("variant_sequence").str().to_uppercase().eq(lit(variant_sequence.to_uppercase())))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_10() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-003-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-003-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-003-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-003-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("position_2").eq(lit(position_2)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_11() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-004-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-004-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-004-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 2);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-004-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 1000)))
            .filter(col("position_1").lt_eq(lit(position_1 + 1000)))
            .filter(col("position_2").gt_eq(lit(position_2 - 1000)))
            .filter(col("position_2").lt_eq(lit(position_2 + 1000)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_12() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-005-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-005-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-005-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-005-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_13() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-006-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-006-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-006-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17","chr18"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-006-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_variant_calling_dna_14() {
    let tumor_bam_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam");
    let tumor_bam_bai_path = Path::new("src/tests/data/bam/dna-007-tumor_minimap2_mdtagged_sorted.bam.bai");
    let normal_bam_path = Path::new("src/tests/data/bam/dna-007-normal_minimap2_mdtagged_sorted.bam");
    let normal_bam_bai_path = Path::new("src/tests/data/bam/dna-007-normal_minimap2_mdtagged_sorted.bam.bai");
    let tumor_bam_full_path = fs::canonicalize(tumor_bam_path).unwrap();
    let tumor_bam_bai_full_path = fs::canonicalize(tumor_bam_bai_path).unwrap();
    let normal_bam_full_path = fs::canonicalize(normal_bam_path).unwrap();
    let normal_bam_bai_full_path = fs::canonicalize(normal_bam_bai_path).unwrap();
    let tumor_bam_file: &str = tumor_bam_full_path.to_str().unwrap();
    let tumor_bam_bai_file: &str = tumor_bam_bai_full_path.to_str().unwrap();
    let normal_bam_file: &str = normal_bam_full_path.to_str().unwrap();
    let normal_bam_bai_file: &str = normal_bam_bai_full_path.to_str().unwrap();
    let control_bam_files: Vec<&str> = vec![normal_bam_file];
    let control_bam_bai_files: Vec<&str> = vec![normal_bam_bai_file];
    let mut variant_call_set: DNAVariantCallSet = identify_case_specific_dna_variants(
        tumor_bam_file,
        tumor_bam_bai_file,
        control_bam_files,
        control_bam_bai_files,
        3,
        30,
        25,
        0.5f32,
        0.5f32,
        2000,
        1000,
        1000,
        true,
        1,
        vec!["chr17","chr18"],
        ""
    );
    assert!(variant_call_set.get_size() == 1);

    // Compare against the ground truth
    let df_variants: DataFrame = variant_call_set.to_dataframe(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/dna-007-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();
        let df_variants_: PolarsResult<DataFrame>;
        df_variants_ = df_variants
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .collect();
        assert!(df_variants_.unwrap().height() == 1);
    }
}

