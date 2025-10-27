use crate::prelude::*;


#[test]
fn test_intron_1() {
    let intron: Intron = Intron::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        1100,
        Strand::Forward,
        1
    );

    assert_eq!(intron.get_size(), 100);
}

#[test]
fn test_intron_2() {
    let intron_1: Intron = Intron::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        1100,
        Strand::Forward,
        1
    );

    let intron_2: Intron = intron_1.clone();

    assert_eq!(intron_1, intron_2);
}
