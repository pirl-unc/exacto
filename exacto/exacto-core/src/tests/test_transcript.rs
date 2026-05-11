use std::collections::HashSet;

use crate::prelude::*;


#[test]
fn test_transcript_1() {
    let mut transcript: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

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

    let exon_2: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE002",
        "havana",
        "chr1",
        1501,
        1600,
        Strand::Forward,
        1,
        2
    );

    let exon_3: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE003",
        "havana",
        "chr1",
        1901,
        2000,
        Strand::Forward,
        1,
        3
    );

    transcript.add_exon(exon_1);
    transcript.add_exon(exon_2);
    transcript.add_exon(exon_3);

    let vectorized_exons = transcript.vectorize_exons(
        "chr1".into(),
        1001,
        2000,
        1,
        0
    );

    assert_eq!(vectorized_exons.len(), 1000);
    assert_eq!(vectorized_exons[0], 1);
    assert_eq!(vectorized_exons[100], 0);
    assert_eq!(vectorized_exons[999], 1);
}

#[test]
fn test_transcript_2() {
    let mut transcript_1: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

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

    let exon_2: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE002",
        "havana",
        "chr1",
        1501,
        1600,
        Strand::Forward,
        1,
        2
    );

    let exon_3: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE003",
        "havana",
        "chr1",
        1901,
        2000,
        Strand::Forward,
        1,
        3
    );

    transcript_1.add_exon(exon_1);
    transcript_1.add_exon(exon_2);
    transcript_1.add_exon(exon_3);

    let transcript_2: Transcript = transcript_1.clone();

    assert_eq!(transcript_1, transcript_2);
    assert_eq!(transcript_2.get_size(), 1000);
    assert_eq!(transcript_2.get_exon("ENSE003").unwrap().exon_number, 3);
    assert_eq!(transcript_2.get_introns().len(), 2);
    assert_eq!(transcript_2.get_introns()[0].start, 1101);
    assert_eq!(transcript_2.get_introns()[0].end, 1500);
    assert_eq!(transcript_2.get_introns()[1].start, 1601);
    assert_eq!(transcript_2.get_introns()[1].end, 1900);
}

#[test]
fn test_transcript_3() {
    let mut transcript: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

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

    let exon_2: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE002",
        "havana",
        "chr1",
        1501,
        1600,
        Strand::Forward,
        1,
        2
    );

    let exon_3: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE003",
        "havana",
        "chr1",
        1901,
        2000,
        Strand::Forward,
        1,
        3
    );

    let utr_5p: UTR = UTR::new(
        "ENSG001",
        "ENST001",
        "ENSE001",
        "havana",
        "chr1",
        1001,
        1020,
        Strand::Forward,
        1,
        1
    );

    let utr_3p: UTR = UTR::new(
        "ENSG001",
        "ENST001",
        "ENSE003",
        "havana",
        "chr1",
        1981,
        2000,
        Strand::Forward,
        1,
        3
    );

    transcript.add_exon(exon_1);
    transcript.add_exon(exon_2);
    transcript.add_exon(exon_3);
    transcript.add_utr(utr_5p);
    transcript.add_utr(utr_3p);

    assert_eq!(transcript.get_5prime_utr().start, 1001);
    assert_eq!(transcript.get_5prime_utr().end, 1020);
    assert_eq!(transcript.get_3prime_utr().start, 1981);
    assert_eq!(transcript.get_3prime_utr().end, 2000);
}

#[test]
fn test_transcript_4() {
    let mut transcript: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

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

    let exon_2: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE002",
        "havana",
        "chr1",
        1501,
        1600,
        Strand::Forward,
        1,
        2
    );

    let exon_3: Exon = Exon::new(
        "ENSG001",
        "ENST001",
        "ENSE003",
        "havana",
        "chr1",
        1901,
        2000,
        Strand::Forward,
        1,
        3
    );

    transcript.add_exon(exon_1);
    transcript.add_exon(exon_2);
    transcript.add_exon(exon_3);

    assert_eq!(transcript.locate_position(1010).0, GenicRegion::Exonic);
    assert_eq!(transcript.locate_position(1010).1.is_some(), true);
    assert_eq!(&*transcript.locate_position(1010).1.unwrap().0.exon_id, "ENSE001");
    assert_eq!(&*transcript.locate_position(1010).1.unwrap().1.exon_id, "ENSE001");

    assert_eq!(transcript.locate_position(1300).0, GenicRegion::Intronic);
    assert_eq!(transcript.locate_position(1300).1.is_some(), true);
    assert_eq!(&*transcript.locate_position(1300).1.unwrap().0.exon_id, "ENSE001");
    assert_eq!(&*transcript.locate_position(1300).1.unwrap().1.exon_id, "ENSE002");

    assert_eq!(transcript.locate_position(3000).0, GenicRegion::Intergenic);
    assert_eq!(transcript.locate_position(3000).1.is_some(), false);
}
