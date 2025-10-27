use crate::prelude::*;


#[test]
fn test_stop_codon_1() {
    let stop_codon: StopCodon = StopCodon::new(
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

    assert_eq!(stop_codon.get_size(), 100);
}

#[test]
fn test_stop_codon_2() {
    let stop_codon_1: StopCodon = StopCodon::new(
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

    let stop_codon_2: StopCodon = stop_codon_1.clone();

    assert_eq!(stop_codon_1, stop_codon_2);
}
