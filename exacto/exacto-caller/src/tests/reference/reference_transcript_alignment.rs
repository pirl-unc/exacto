use bio::alignment::AlignmentMode;
use exacto_core::prelude::Nucleotide;
use std::str::FromStr;

use crate::prelude::*;

use super::*;


/// A 60-base transcript, arbitrary but non-repetitive.
const TRANSCRIPT: &str = "ACGTTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAGTTCAGCATGGACTTGCACGAAT";


fn build_reference_transcript_sequence(sequence: &str) -> ReferenceTranscriptSequence {
    let mut rts: ReferenceTranscriptSequence = ReferenceTranscriptSequence::new(
        "GENE",
        "TRANSCRIPT"
    );

    for (i, base) in sequence.chars().enumerate() {
        rts.push(ReferenceBase::new(
            0,
            1_000 + i as u32,
            Nucleotide::from_str(&base.to_string()).unwrap(),
            Strand::Forward,
            None,
            None,
            None
        ));
    }

    rts
}


fn place(
    query: &str,
    rts: &ReferenceTranscriptSequence,
    window: Range<u32>,
    min_score_fraction: f64,
    min_query_coverage: f64
) -> Option<ReferenceTranscriptPlacement> {
    place_on_reference_transcript(
        query,
        rts,
        window,
        -5,
        -1,
        5,
        10,
        min_score_fraction,
        min_query_coverage
    )
}

/// Complement the bases at `positions`, which guarantees a mismatch at each one.
fn complement_at(sequence: &str, positions: &[usize]) -> String {
    sequence
        .chars()
        .enumerate()
        .map(|(i, base)| if positions.contains(&i) {
            match base {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                other => other
            }
        } else {
            base
        })
        .collect()
}


#[test]
fn reference_transcript_alignment_places_forward_exact_match() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let placement: ReferenceTranscriptPlacement = place(
        &TRANSCRIPT[10..40],
        &rts,
        0..60,
        0.7,
        0.5
    ).unwrap();

    assert!(placement.is_forward);
    assert_eq!(placement.query_range, 0..30);
    assert_eq!(placement.reference_range, 10..40);
    assert_eq!(placement.score, 30);
}


#[test]
fn reference_transcript_alignment_places_reverse_match() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = reverse_complement(&TRANSCRIPT[10..40]).to_string();
    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        0.7,
        0.5
    ).unwrap();

    assert!(!placement.is_forward);
    assert_eq!(placement.query_range, 0..30);
    assert_eq!(placement.reference_range, 10..40);
    assert_eq!(placement.score, 30);
}


#[test]
fn reference_transcript_alignment_offsets_reference_range_by_window_start() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let placement: ReferenceTranscriptPlacement = place(
        &TRANSCRIPT[25..50],
        &rts,
        20..60,
        0.7,
        0.5
    ).unwrap();

    assert!(placement.is_forward);
    assert_eq!(placement.query_range, 0..25);
    assert_eq!(placement.reference_range, 25..50);
    assert_eq!(placement.score, 25);
}


#[test]
fn reference_transcript_alignment_rejects_query_outside_window() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = reverse_complement(&TRANSCRIPT[35..60]).to_string();
    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..30,
        0.7,
        0.5
    );

    assert!(placement.is_none());
}


/// 12 aligned bases of a 28-base query is 0.43 coverage: rejected at 0.5, accepted at 0.4.
#[test]
fn reference_transcript_alignment_coverage_gate() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = format!("{}{}", "N".repeat(16), &TRANSCRIPT[10..22]).to_string();
    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..60,
        0.7,
        0.5
    );

    assert!(placement.is_none());

    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        0.7,
        0.4
    ).unwrap();

    assert!(placement.is_forward);
    assert_eq!(placement.query_range, 16..28);
    assert_eq!(placement.reference_range, 10..22);
    assert_eq!(placement.score, 12);
}


/// Three double mismatches inside a 26-base query: 20 matches, 6 mismatches, score 14,
/// fraction 0.54. Rejected at 0.7, accepted at 0.5. The runs between mismatches are five
/// bases so every block still carries a 5-mer seed; single mismatches every fifth base
/// would leave the banded aligner nothing to seed on.
#[test]
fn reference_transcript_alignment_score_fraction_gate() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = complement_at(&TRANSCRIPT[10..36], &[5, 6, 12, 13, 19, 20]);
    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..60,
        0.7,
        0.5
    );

    assert!(placement.is_none());

    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        0.5,
        0.5
    ).unwrap();

    assert!(placement.is_forward);
    assert_eq!(placement.query_range, 0..26);
    assert_eq!(placement.reference_range, 10..36);
    assert_eq!(placement.score, 14);
}


/// Both strands clear the gates; the higher score wins.
/// Forward matches 10 bases while the reverse complement matches 30.
#[test]
fn reference_transcript_alignment_prefers_higher_scoring_strand() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = format!("{}{}", &TRANSCRIPT[0..10], reverse_complement(&TRANSCRIPT[20..50]));
    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        0.7,
        0.2
    ).unwrap();

    assert!(!placement.is_forward);
    assert_eq!(placement.query_range, 0..30);
    assert_eq!(placement.reference_range, 20..50);
    assert_eq!(placement.score, 30);
}


/// A palindromic query scores the same in both orientations. The comparison is strict, so
/// the forward candidate, tried first, is kept.
#[test]
fn reference_transcript_alignment_equal_orientation_scores_keep_forward() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence("ATGCAT");
    let placement: ReferenceTranscriptPlacement = place(
        "ATGCAT",
        &rts,
        0..6,
        1.0,
        1.0
    ).unwrap();

    assert!(placement.is_forward);
    assert_eq!(placement.query_range, 0..6);
    assert_eq!(placement.reference_range, 0..6);
    assert_eq!(placement.score, 6);
}


/// Eight aligned bases of a ten-base query is exactly 0.8 coverage: accepted at 0.8,
/// rejected one hundredth above it.
#[test]
fn reference_transcript_alignment_coverage_threshold_is_inclusive() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = format!("NN{}", &TRANSCRIPT[10..18]);
    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        1.0,
        0.8
    ).unwrap();

    assert_eq!(placement.query_range, 2..10);
    assert_eq!(placement.reference_range, 10..18);
    assert_eq!(placement.score, 8);

    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..60,
        1.0,
        0.81
    );

    assert!(placement.is_none());
}


/// Seven matches and one mismatch over eight aligned bases is exactly 0.75: accepted at
/// 0.75, rejected one hundredth above it.
#[test]
fn reference_transcript_alignment_score_fraction_threshold_is_inclusive() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let query: String = complement_at(&TRANSCRIPT[10..18], &[4]);
    let placement: ReferenceTranscriptPlacement = place(
        &query,
        &rts,
        0..60,
        0.75,
        1.0
    ).unwrap();

    assert_eq!(placement.query_range, 0..8);
    assert_eq!(placement.reference_range, 10..18);
    assert_eq!(placement.score, 6);

    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..60,
        0.76,
        1.0
    );

    assert!(placement.is_none());
}


/// Neither CCCCCCCC nor its reverse complement shares a base with the transcript, so the
/// alignment has zero length. That is rejected before the gates, which are both zero here.
#[test]
fn reference_transcript_alignment_no_alignable_bases_returns_none_with_zero_thresholds() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence("AAAAAAAA");
    let query: String = "CCCCCCCC".to_string();
    let placement: Option<ReferenceTranscriptPlacement> = place(
        &query,
        &rts,
        0..8,
        0.0,
        0.0
    );

    assert!(placement.is_none());
}


/// An empty query already returns `None`. An empty reference window never returned at all
/// until the guard at the top of `place_on_reference_transcript` was added: rust-bio's
/// banded traceback loops forever when the reference is empty.
#[test]
fn reference_transcript_alignment_empty_inputs_return_none() {
    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence(TRANSCRIPT);
    let placement: Option<ReferenceTranscriptPlacement> = place(
        "",
        &rts,
        0..60,
        0.0,
        0.0
    );

    assert!(placement.is_none());

    let placement: Option<ReferenceTranscriptPlacement> = place(
        "ACGT",
        &rts,
        4..4,
        0.0,
        0.0
    );

    assert!(placement.is_none());

    let rts: ReferenceTranscriptSequence = build_reference_transcript_sequence("");
    let placement: Option<ReferenceTranscriptPlacement> = place(
        "ACGT",
        &rts,
        0..0,
        0.0,
        0.0
    );

    assert!(placement.is_none());
}
#[test]
fn reference_transcript_alignment_places_a_query_on_a_soft_masked_transcript() {
    // The FASTA holds the transcript in lowercase (soft-masked); read bases are uppercase.
    let transcript: &str = "acgttgcaagcttaggcatccgataagtcgtacctgagttcagcatggacttgcacgaat";
    let query: String = transcript[10..40].to_ascii_uppercase();
    let mut rts: ReferenceTranscriptSequence = ReferenceTranscriptSequence::new("GENE", "TRANSCRIPT");
    for (i, base) in transcript.chars().enumerate() {
        rts.push(ReferenceBase::new(
            0,
            1_000 + i as u32,
            Nucleotide::from_str(&base.to_string()).unwrap(),
            Strand::Forward,
            None,
            None,
            None
        ));
    }

    let placement: ReferenceTranscriptPlacement = place_on_reference_transcript(&query, &rts, 0..60, -5, -1, 5, 10, 0.8, 0.5)
        .expect("the query is placed on the transcript whatever the case of its bases");

    assert_eq!(placement.query_range, 0..30);
    assert_eq!(placement.reference_range, 10..40);
    assert!(placement.is_forward);
    assert_eq!(placement.score, 30);
}
