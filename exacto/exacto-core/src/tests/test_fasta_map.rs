use std::fs;
use std::path::Path;

use crate::index::fasta_map::FastaMap;


#[test]
fn test_fasta_map() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);
    let sequence: &str = fasta_map.get_sequence("chr17", 7668784, 7668789);
    assert!(sequence == "GCCACT");
}

#[test]
fn test_fasta_map_soft_masked() {
    let temp_dir = tempfile::TempDir::new().unwrap();
    let fasta_file: String = temp_dir.path().join("masked.fa").to_str().unwrap().to_string();
    fs::write(&fasta_file, ">chrM\nACGTacgtNNnnACGT\n").unwrap();
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    assert_eq!(fasta_map.get_sequence("chrM", 1, 16), "ACGTACGTNNNNACGT");
}
