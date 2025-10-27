use crate::prelude::*;

#[test]
fn test_cds_1() {
    let cds: CDS = CDS::new(
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
    assert_eq!(cds.get_size(), 100);
}

#[test]
fn test_cds_2() {
    let cds_1: CDS = CDS::new(
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

    let cds_2: CDS = cds_1.clone();

    assert_eq!(cds_1, cds_2);
}
