use exacto_core::prelude::*;
use noodles_bam as bam;
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::prelude::*;


#[test]
fn test_remove_unspliced_rnas_1() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bai_full_path = fs::canonicalize(bai_path).unwrap();
    let bai_file: &str = bai_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    remove_unspliced_rnas(
        bam_file,
        bai_file,
        &gencode,
        output_bam_file.as_str(),
        output_bai_file.as_str(),
        1,
        30,
        true
    );

    // 358 of the 403 reads keep a primary or supplementary record at MAPQ >= 30 that splices or
    // overlaps a single-exon transcript; each keeps one record.
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(output_bam_file.as_str())
        .unwrap();
    reader.read_header().unwrap();
    let num_records: usize = reader.records().filter_map(Result::ok).count();
    assert_eq!(num_records, 358);
}

#[test]
fn test_remove_unspliced_rnas_decides_each_read_on_one_record() {
    // chrS (nothing annotated there, so no read is kept by a single-exon transcript):
    //   ctrl_spliced                primary 50M200N50M MAPQ 60                         -> kept
    //   ctrl_unspliced              primary 100M MAPQ 60                               -> removed
    //   secondary_spliced           primary 100M MAPQ 60, secondary 50M200N50M MAPQ 60 -> removed
    //   split_criteria              primary 50M200N50M50S MAPQ 0,
    //                               supplementary 100H50M MAPQ 60                      -> removed
    //   spliced_with_other_records  primary 50M200N50M50S MAPQ 60,
    //                               supplementary 100H50M MAPQ 60, secondary 150M     -> kept,
    //                               without its secondary record
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-qc/unspliced_rna_filtering_records.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-qc/unspliced_rna_filtering_records.bam.bai");
    let bai_full_path = fs::canonicalize(bai_path).unwrap();
    let bai_file: &str = bai_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let read_names_to_keep: HashSet<Box<str>> = remove_unspliced_rnas(
        bam_file,
        bai_file,
        &gencode,
        output_bam_file.as_str(),
        output_bai_file.as_str(),
        1,
        30,
        true
    );

    assert_eq!(read_names_to_keep, HashSet::from(["ctrl_spliced".into(), "spliced_with_other_records".into()]));
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(output_bam_file.as_str())
        .unwrap();
    reader.read_header().unwrap();
    let written: Vec<(String, u16, usize)> = reader
        .records()
        .map(|result| {
            let record: bam::Record = result.unwrap();
            (
                record.name().unwrap().to_string(),
                record.flags().bits(),
                record.alignment_start().unwrap().unwrap().get()
            )
        })
        .collect();
    assert_eq!(written, vec![
        ("ctrl_spliced".to_string(), 0, 101),
        ("spliced_with_other_records".to_string(), 0, 1501),
        ("spliced_with_other_records".to_string(), 2048, 4201),
    ]);
}


#[test]
fn test_remove_unspliced_rnas_refuses_to_overwrite_its_input() {
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    // The input is a copy, named once plainly and once through a `.` component.
    let temp_dir = TempDir::new().unwrap();
    let bam_file: String = temp_dir.path().join("input.bam").to_str().unwrap().to_string();
    let bai_file: String = temp_dir.path().join("input.bam.bai").to_str().unwrap().to_string();
    let same_bam_file: String = temp_dir.path().join(".").join("input.bam").to_str().unwrap().to_string();
    let other_bai_file: String = temp_dir.path().join("output.bam.bai").to_str().unwrap().to_string();
    fs::copy(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-qc/unspliced_rna_filtering_records.bam"), &bam_file).unwrap();
    fs::copy(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-qc/unspliced_rna_filtering_records.bam.bai"), &bai_file).unwrap();
    let bam_bytes: Vec<u8> = fs::read(&bam_file).unwrap();
    let bai_bytes: Vec<u8> = fs::read(&bai_file).unwrap();

    for (output_bam_file, output_bai_file) in [
        (same_bam_file.as_str(), other_bai_file.as_str()),
        (temp_dir.path().join("output.bam").to_str().unwrap(), bai_file.as_str()),
    ] {
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            remove_unspliced_rnas(
                bam_file.as_str(),
                bai_file.as_str(),
                &gencode,
                output_bam_file,
                output_bai_file,
                1,
                30,
                true
            )
        }));
        assert!(result.is_err(), "writing {output_bam_file} and {output_bai_file} was not refused");
        assert_eq!(fs::read(&bam_file).unwrap(), bam_bytes);
        assert_eq!(fs::read(&bai_file).unwrap(), bai_bytes);
    }
}
