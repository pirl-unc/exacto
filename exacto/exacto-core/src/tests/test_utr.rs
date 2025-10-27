use crate::prelude::*;


#[test]
fn test_utr_1() {
    let utr: UTR = UTR::new(
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

    assert_eq!(utr.get_size(), 100);
}

#[test]
fn test_utr_2() {
    let utr_1: UTR = UTR::new(
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

    let utr_2: UTR = utr_1.clone();

    assert_eq!(utr_1, utr_2);
}
