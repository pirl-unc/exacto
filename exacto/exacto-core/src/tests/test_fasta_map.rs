use std::fs;
use std::path::Path;

use crate::index::fasta_map::FastaMap;


#[test]
fn test_fasta_map() {
    let fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);
    let sequence: &str = fasta_map.get_sequence("chr17", 7668784, 7668789);
    assert!(sequence == "gccact");
}
