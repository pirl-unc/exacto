use crate::prelude::*;
use crate::common::constants::Strands;


#[test]
fn test_vectorize_exons() {
    let mut transcript: Transcript = Transcript::new(
        "gene_1",
        "transcript_1",
        "custom",
        "chr1",
        3501,
        3600,
        Strands::Forward,
        1,
        "transcript_1",
        "custom",
        "1"
    );

    let exon_1: Exon = Exon::new(
        "gene_1",
        "transcript_1",
        "exon_1",
        "custom",
        "chr1",
        3501,
        3510,
        Strands::Forward,
        1,
        1
    );

    let exon_2: Exon = Exon::new(
        "gene_1",
        "transcript_1",
        "exon_2",
        "custom",
        "chr1",
        3551,
        3560,
        Strands::Forward,
        1,
        2
    );

    let exon_3: Exon = Exon::new(
        "gene_1",
        "transcript_1",
        "exon_3",
        "custom",
        "chr1",
        3591,
        3600,
        Strands::Forward,
        1,
        3
    );

    transcript.add_exon(exon_1);
    transcript.add_exon(exon_2);
    transcript.add_exon(exon_3);

    let vectorized_exons = transcript.vectorize_exons(
        "chr1".into(),
        3501,
        3600,
        1,
        0
    );

    assert_eq!(vectorized_exons.len(), 100);
    assert_eq!(vectorized_exons[0], 1);
    assert_eq!(vectorized_exons[20], 0);
    assert_eq!(vectorized_exons[99], 1);
}