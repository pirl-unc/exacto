use crate::prelude::*;


#[test]
fn test_start_codon_1() {
    let start_codon: StartCodon = StartCodon::new(
        "ENSG001",
        "ENST001",
        "ENSE001",
        "havana",
        "chr1",
        1001,
        1100,
        Strand::Forward,
        1,
        1
    );

    assert_eq!(start_codon.get_size(), 100);
}

#[test]
fn test_start_codon_2() {
    let start_codon_1: StartCodon = StartCodon::new(
        "ENSG001",
        "ENST001",
        "ENSE001",
        "havana",
        "chr1",
        1001,
        1100,
        Strand::Forward,
        1,
        1
    );

    let start_codon_2: StartCodon = start_codon_1.clone();

    assert_eq!(start_codon_1, start_codon_2);
}
