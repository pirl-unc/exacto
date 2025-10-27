use bio::io::bed;
use polars::prelude::*;
use std::fs;
use std::path::Path;

use crate::common::files::*;


#[test]
fn test_read_tsv_file_1() {
    let tsv_path = Path::new("src/tests/data/tsv/samples.tsv.gz");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();
    let df: DataFrame = read_tsv_file(tsv_file);
    assert!(df.height() == 1);
}

#[test]
fn test_read_fasta_file_1() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let sequences: Vec<(Box<str>,Box<str>)> = read_fasta_file(fasta_file);
    assert!(sequences.len() == 2);
}

#[test]
fn test_read_fasta_file_2() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let sequences: Vec<(Box<str>,Box<str>)> = read_fasta_file(fasta_file);
    assert!(sequences.len() == 2);
}

#[test]
fn test_read_bed_file_1() {
    let path = Path::new("src/tests/data/bed/hg38_chr17-18.bed");
    let full_path = fs::canonicalize(path).unwrap();
    let bed_file: &str = full_path.to_str().unwrap();
    let records: Vec<bed::Record> = read_bed_file(bed_file);
    assert!(records.len() == 2);
}

