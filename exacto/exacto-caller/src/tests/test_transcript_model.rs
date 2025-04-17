use crate::prelude::*;


#[test]
fn test_transcript_model_vectorize_exons_1() {
    let mut transcript_model = TranscriptModel::new(
        0,
        Vec::new(),
        Vec::new()
    );
    
    let exon_1: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3501,
        3510,
        1,
        Strands::Forward
    );

    let exon_2: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3551,
        3560,
        2,
        Strands::Forward
    );

    let exon_3: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3591,
        3600,
        3,
        Strands::Forward
    );

    transcript_model.add_exon(exon_1);
    transcript_model.add_exon(exon_2);
    transcript_model.add_exon(exon_3);
    
    let vectorized_exons = transcript_model.vectorize_exons(
        1,
        3501,
        3600,
        1,
        0
    );
    
    assert!(vectorized_exons.len() == 100);
    assert!(vectorized_exons[0] == 1);
    assert!(vectorized_exons[20] == 0);
    assert!(vectorized_exons[99] == 1);
}

#[test]
fn test_transcript_model_vectorize_exons_2() {
    let mut transcript_model = TranscriptModel::new(
        0,
        Vec::new(),
        Vec::new()
    );

    let exon_1: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3501,
        3510,
        1,
        Strands::Forward
    );

    let exon_2: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3551,
        3560,
        2,
        Strands::Forward
    );

    let exon_3: TranscriptModelExon = TranscriptModelExon::new(
        1,
        3591,
        3600,
        3,
        Strands::Forward
    );

    transcript_model.add_exon(exon_1);
    transcript_model.add_exon(exon_2);
    transcript_model.add_exon(exon_3);

    let vectorized_exons = transcript_model.vectorize_exons(
        1,
        3599,
        3601,
        1,
        0
    );

    assert!(vectorized_exons.len() == 3);
    assert!(vectorized_exons[0] == 1);
    assert!(vectorized_exons[1] == 1);
    assert!(vectorized_exons[2] == 0);
}