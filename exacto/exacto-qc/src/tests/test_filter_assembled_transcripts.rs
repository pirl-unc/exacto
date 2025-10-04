use exacto_core::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::prelude::*;


#[test]
fn test_filter_assembled_transcripts() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted_dna-001-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bai_full_path = fs::canonicalize(bai_path).unwrap();
    let bai_file: &str = bai_full_path.to_str().unwrap();
    let fasta_path = Path::new("src/tests/data/fasta/rna-100-tumor_long-read_dna-001-tumor_long-read.fasta");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let gencode: Gencode = Gencode::new(
        gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from_iter(vec!["protein_coding"])),
        Some(HashSet::from_iter(vec![1,2])),
        Some(HashSet::from_iter(vec!["protein_coding"])),
        Some(HashSet::from_iter(vec![1,2]))
    );

    let temp_dir = TempDir::new().unwrap();
    let output_bam_file: String = temp_dir.path().join("test.bam").to_str().unwrap().to_string();
    let output_bai_file: String = temp_dir.path().join("test.bam.bai").to_str().unwrap().to_string();
    let output_fasta_file: String = temp_dir.path().join("test.fasta").to_str().unwrap().to_string();

    filter_assembled_transcripts(
        bam_file,
        bai_file,
        fasta_file,
        &gencode,
        output_bam_file.as_str(),
        output_bai_file.as_str(),
        output_fasta_file.as_str(),
        1,
        30
    );

    assert!(get_fasta_sequence_ids(output_fasta_file.as_str()).len() == 2);
}
