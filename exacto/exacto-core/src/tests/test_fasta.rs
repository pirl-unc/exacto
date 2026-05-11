use polars::prelude::*;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::common::fasta::*;


#[test]
fn test_create_fai_file() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    create_fai_file(fasta_file);
}

#[test]
fn test_fasta_index_exists() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    assert!(fasta_index_exists(fasta_file) == true);
}

#[test]
fn test_get_fasta_sequence() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let sequence: Box<str> = get_fasta_sequence("chr17", 7668784, 7668798, fasta_file);
    assert!(&*sequence == "gccactccactccag");
}

#[test]
fn test_get_fasta_sequence_ids() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let sequence_ids: Vec<(Box<str>, u32)> = get_fasta_sequence_ids(fasta_file);
    assert!(sequence_ids.len() == 2);
    assert!(sequence_ids[0].0 == "chr17".into());
    assert!(sequence_ids[1].0 == "chr18".into());
}

#[test]
fn test_write_fasta_file() {
    let temp_dir = TempDir::new().unwrap();
    let fasta_file: String = temp_dir.path().join("test.fasta").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("seq1".into(), "TCGA".into()), 
        ("seq2".into(), "ACTG".into())
    ];
    write_fasta_file(&sequences, fasta_file.as_str());
}

