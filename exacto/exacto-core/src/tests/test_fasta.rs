use polars::prelude::*;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::common::fasta::*;


#[test]
fn test_create_fai_file() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    create_fai_file(fasta_file);
}

#[test]
fn test_fasta_index_exists() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    assert!(fasta_index_exists(fasta_file) == true);
}

#[test]
fn test_get_fasta_sequence() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let sequence: Box<str> = get_fasta_sequence("chr17", 7668784, 7668798, fasta_file);
    assert!(&*sequence == "GCCACTCCACTCCAG");
}

#[test]
fn test_get_fasta_sequence_soft_masked() {
    // Soft-masked bases (lowercase, as in UCSC hg38) come back in uppercase, as read bases are.
    let temp_dir = TempDir::new().unwrap();
    let fasta_file: String = temp_dir.path().join("masked.fa").to_str().unwrap().to_string();
    fs::write(&fasta_file, ">chrM\nACGTacgtNNnnACGT\n").unwrap();
    assert_eq!(&*get_fasta_sequence("chrM", 1, 16, &fasta_file), "ACGTACGTNNNNACGT");
    assert_eq!(&*get_fasta_sequence("chrM", 5, 8, &fasta_file), "ACGT");
}

#[test]
fn test_get_fasta_sequence_parallel_first_use() {
    // Threads that open a FASTA without an index at the same moment create the index once and
    // all read their sequence; none reads an index half written.
    let temp_dir = TempDir::new().unwrap();
    for file_number in 0..5 {
        let fasta_file: String = temp_dir.path().join(format!("parallel_{file_number}.fa")).to_str().unwrap().to_string();
        let mut content: String = String::new();
        for contig in 0..300 {
            content.push_str(&format!(">contig{contig}\n{}\n", "ACGT".repeat(20)));
        }
        fs::write(&fasta_file, content).unwrap();

        let barrier = std::sync::Arc::new(std::sync::Barrier::new(16));
        let handles: Vec<std::thread::JoinHandle<Box<str>>> = (0..16)
            .map(|_| {
                let barrier = barrier.clone();
                let fasta_file: String = fasta_file.clone();
                std::thread::spawn(move || {
                    barrier.wait();
                    get_fasta_sequence("contig299", 1, 8, &fasta_file)
                })
            })
            .collect();
        for handle in handles {
            assert_eq!(&*handle.join().unwrap(), "ACGTACGT");
        }
    }
}

#[test]
fn test_get_fasta_sequence_ids() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
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

