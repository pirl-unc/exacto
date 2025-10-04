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

