use exacto_core::prelude::*;

use crate::prelude::*;


#[test]
fn test_assembled_transcript_exon_1() {
    let exon_1: AssembledTranscriptExon = AssembledTranscriptExon::new(
        0,
        1,
        30,
        Strand::Forward,
        1,
        0,
        29
    );
    let exon_2: AssembledTranscriptExon = AssembledTranscriptExon::new(
        0,
        71,
        100,
        Strand::Forward,
        2,
        30,
        59
    );
    let mut exons: Vec<AssembledTranscriptExon> = Vec::new();
    exons.push(exon_1);
    exons.push(exon_2);
    
    let vectorized_exons: Vec<i8> = vectorize_exons(
        &exons,
        0,
        1,
        100,
        1i8,
        0i8
    );

    assert!(vectorized_exons.len() == 100);
    assert!(vectorized_exons[0] == 1);
    assert!(vectorized_exons[31] == 0);
    assert!(vectorized_exons[99] == 1);
}
