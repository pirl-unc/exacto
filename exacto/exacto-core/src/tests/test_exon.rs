use crate::prelude::*;

#[test]
fn test_exon_1() {
    let exon: Exon = Exon::new(
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

    assert_eq!(exon.get_size(), 100);
}

#[test]
fn test_exon_2() {
    let exon_1: Exon = Exon::new(
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

    let exon_2: Exon = exon_1.clone();

    assert_eq!(exon_1, exon_2);
}
