use exacto_core::prelude::{reverse_complement, Nucleotide};
use std::str::FromStr;

use crate::prelude::*;

use super::*;


/// A 60-base transcript, arbitrary but non-repetitive.
const TRANSCRIPT: &str = "ACGTTGCAAGCTTAGGCATCCGATAAGTCGTACCTGAGTTCAGCATGGACTTGCACGAAT";

const CHROMOSOME: u16 = 0;
const READ_ID: usize = 7;
const READ_POSITION_1: u32 = 100;
const READ_POSITION_2: u32 = 101;


/// Transcript base `i` sits at the `i`-th of `positions`. Forward transcripts get ascending
/// positions and reverse transcripts descending ones, matching how
/// `ReferenceTranscriptSequence::from_reference_transcript` lays out exons.
fn build_transcript_helper(
    chromosome: u16,
    sequence: &str,
    positions: impl IntoIterator<Item = u32>,
    strand: Strand
) -> ReferenceTranscriptSequence {
    let mut rts: ReferenceTranscriptSequence = ReferenceTranscriptSequence::new(
        "GENE",
        "TRANSCRIPT"
    );

    for (base, position) in sequence.chars().zip(positions) {
        rts.push(ReferenceBase::new(
            chromosome,
            position,
            Nucleotide::from_str(&base.to_string()).unwrap(),
            strand.clone(),
            None,
            None,
            None
        ));
    }

    assert_eq!(rts.get_length(), sequence.len(), "every base needs a position");
    rts
}


/// A transcript on the default chromosome.
fn build_transcript(
    sequence: &str,
    positions: impl IntoIterator<Item = u32>,
    strand: Strand
) -> ReferenceTranscriptSequence {
    build_transcript_helper(CHROMOSOME, sequence, positions, strand)
}


/// An insertion of `sequence` between the two flanking genomic positions.
fn insertion(position_1: u32, position_2: u32, strand: Strand, sequence: &str) -> VariantRecord {
    VariantRecord::new(
        READ_ID,
        READ_POSITION_1,
        READ_POSITION_2,
        GraphOperation::new(
            CHROMOSOME,
            position_1,
            strand.clone(),
            GraphOperationType::Include,
            CHROMOSOME,
            position_2,
            strand,
            GraphOperationType::Include,
            sequence.into(),
            VariantType::Insertion
        )
    )
}


/// The breakpoint record the rescue emits. `GraphOperation::new` standardizes to ascending
/// position, so every expected record below is written in that order.
fn breakpoint(
    position_1: u32,
    strand_1: Strand,
    operation_1: GraphOperationType,
    position_2: u32,
    strand_2: Strand,
    operation_2: GraphOperationType
) -> VariantRecord {
    VariantRecord::new(
        READ_ID,
        READ_POSITION_1,
        READ_POSITION_2,
        GraphOperation::new(
            CHROMOSOME,
            position_1,
            strand_1,
            operation_1,
            CHROMOSOME,
            position_2,
            strand_2,
            operation_2,
            "".into(),
            VariantType::Breakpoint
        )
    )
}


fn retype(
    record: &VariantRecord,
    transcripts: &Vec<ReferenceTranscriptSequence>,
    min_placed_fraction: f64,
    max_pieces: u32
) -> Option<Vec<VariantRecord>> {
    retype_insertion_on_reference_transcripts(
        record,
        1_000,
        transcripts,
        10,
        -5,
        -1,
        5,
        10,
        0.8,
        0.3,
        min_placed_fraction,
        max_pieces
    )
}


#[test]
#[should_panic]
fn breakpoint_rescue_panics_on_a_record_that_is_not_an_insertion() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let deletion: VariantRecord = VariantRecord::new(
        READ_ID,
        READ_POSITION_1,
        READ_POSITION_2,
        GraphOperation::new(
            CHROMOSOME,
            1_029,
            Strand::Forward,
            GraphOperationType::Downstream,
            CHROMOSOME,
            1_040,
            Strand::Forward,
            GraphOperationType::Upstream,
            "".into(),
            VariantType::Deletion
        )
    );

    retype(&deletion, &transcripts, 0.0, 3);
}


#[test]
fn breakpoint_rescue_length_gate_is_inclusive() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];

    // Nine bases: one short of the ten-base minimum, even though they would place perfectly.
    let short: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..19]);
    assert!(retype(&short, &transcripts, 0.0, 3).is_none());

    // Ten bases clear the gate.
    let minimum: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..20]);
    assert!(retype(&minimum, &transcripts, 0.0, 3).is_some());
}


#[test]
fn breakpoint_rescue_returns_none_when_no_transcript_holds_both_flanks() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];

    // Right flank one base past the transcript end.
    let record: VariantRecord = insertion(1_059, 1_060, Strand::Forward, &TRANSCRIPT[10..30]);
    assert!(retype(&record, &transcripts, 0.0, 3).is_none());

    // Neither flank on the transcript.
    let record: VariantRecord = insertion(2_029, 2_030, Strand::Forward, &TRANSCRIPT[10..30]);
    assert!(retype(&record, &transcripts, 0.0, 3).is_none());
}


/// The insertion copies transcript bases 10..30 and sits right after base 29. The walk
/// leaves base 29 for base 10, a back junction, and then leaves base 29 for base 30, which
/// is ordinary continuation and is dropped. The read fields are carried over unchanged.
#[test]
fn breakpoint_rescue_tandem_duplication_yields_one_back_junction() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// The insertion copies bases 30..50, the bases that follow the insertion point, so the
/// first hop (29 -> 30) is continuation and only the return hop (49 -> 30) is a junction.
#[test]
fn breakpoint_rescue_drops_leading_continuation() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[30..50]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_030, Strand::Forward, GraphOperationType::Upstream,
            1_049, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// A copy of bases 40..55 inserted after base 9 touches neither flank, so both hops are
/// junctions, reported in walk order: 9 -> 40, then 54 -> 10.
#[test]
fn breakpoint_rescue_transposed_copy_yields_two_junctions_in_walk_order() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_009, 1_010, Strand::Forward, &TRANSCRIPT[40..55]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_009, Strand::Forward, GraphOperationType::Downstream,
            1_040, Strand::Forward, GraphOperationType::Upstream
        ),
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_054, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// The reverse complement of bases 10..30 inserted after base 29. The walk enters the
/// piece at its transcript end (base 29) on the flipped strand and leaves it at base 10,
/// still flipped, before rejoining base 30 on the transcript strand.
#[test]
fn breakpoint_rescue_inverted_duplication_flips_strand_across_the_piece() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = reverse_complement(&TRANSCRIPT[10..30]).to_string();
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_029, Strand::Forward, GraphOperationType::Downstream,
            1_029, Strand::Reverse, GraphOperationType::Downstream
        ),
        breakpoint(
            1_010, Strand::Reverse, GraphOperationType::Upstream,
            1_030, Strand::Forward, GraphOperationType::Upstream
        )
    ]);
}


/// Bases 40..52 followed by bases 2..20. The 18-base piece is placed first because it
/// scores higher, but the chain follows read order: 29 -> 40, 51 -> 2, 19 -> 30.
#[test]
fn breakpoint_rescue_chains_pieces_in_read_order_not_placement_order() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = format!("{}{}", &TRANSCRIPT[40..52], &TRANSCRIPT[2..20]);
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_029, Strand::Forward, GraphOperationType::Downstream,
            1_040, Strand::Forward, GraphOperationType::Upstream
        ),
        breakpoint(
            1_002, Strand::Forward, GraphOperationType::Upstream,
            1_051, Strand::Forward, GraphOperationType::Downstream
        ),
        breakpoint(
            1_019, Strand::Forward, GraphOperationType::Downstream,
            1_030, Strand::Forward, GraphOperationType::Upstream
        )
    ]);
}


/// Bases 40..52 forward, then the reverse complement of bases 5..20. The reverse piece is
/// found first, at 0..15 on the reverse-complemented query. That must map back to 12..27
/// on the original query so the mask leaves the forward piece intact and the chain stays
/// in read order: 29 -> 40 forward, 51 -> 19 entering the piece backwards, 5 -> 30.
#[test]
fn breakpoint_rescue_maps_reverse_piece_back_to_query_coordinates() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = format!(
        "{}{}",
        &TRANSCRIPT[40..52],
        reverse_complement(&TRANSCRIPT[5..20])
    );
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_029, Strand::Forward, GraphOperationType::Downstream,
            1_040, Strand::Forward, GraphOperationType::Upstream
        ),
        breakpoint(
            1_019, Strand::Reverse, GraphOperationType::Downstream,
            1_051, Strand::Forward, GraphOperationType::Downstream
        ),
        breakpoint(
            1_005, Strand::Reverse, GraphOperationType::Upstream,
            1_030, Strand::Forward, GraphOperationType::Upstream
        )
    ]);
}


/// Twenty placeable bases followed by ten Ns place 0.67 of the insertion: rejected at 0.9,
/// accepted at 0.6 with the trailing Ns simply left out of the chain.
#[test]
fn breakpoint_rescue_placed_fraction_gate() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = format!("{}{}", &TRANSCRIPT[10..30], "N".repeat(10));
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    assert!(retype(&record, &transcripts, 0.9, 3).is_none());

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.6, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// With one piece allowed, only the higher-scoring 18-base piece is placed: 0.6 of the
/// insertion, so 0.9 rejects it and 0.5 accepts a chain through that piece alone.
#[test]
fn breakpoint_rescue_max_pieces_caps_the_decomposition() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = format!("{}{}", &TRANSCRIPT[40..52], &TRANSCRIPT[2..20]);
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    assert!(retype(&record, &transcripts, 0.9, 1).is_none());

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.5, 1).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_002, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        ),
        breakpoint(
            1_019, Strand::Forward, GraphOperationType::Downstream,
            1_030, Strand::Forward, GraphOperationType::Upstream
        )
    ]);
}


/// With nothing placed and a zero threshold the only hop is 29 -> 30, plain continuation.
/// The rescue must answer `None`, not `Some(vec![])`: the caller replaces the insertion
/// with whatever comes back, and an empty list would silently delete the insertion.
#[test]
fn breakpoint_rescue_never_returns_an_empty_record_list() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];

    // Nothing to place.
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &"N".repeat(20));
    assert!(retype(&record, &transcripts, 0.0, 3).is_none());

    // Placeable, but no pieces allowed.
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);
    assert!(retype(&record, &transcripts, 0.0, 0).is_none());
}


/// The same tandem duplication on a reverse-strand transcript, whose bases run from 2_059
/// down to 2_000. Transcript offsets are unchanged, but the junction is reported on the
/// reverse strand with the orientations swapped: Upstream at 2_030, Downstream at 2_049.
#[test]
fn breakpoint_rescue_reports_reverse_strand_transcript_junctions_on_the_reverse_strand() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, (2_000..2_060).rev(), Strand::Reverse)
    ];

    // Flanks are transcript offsets 29 and 30, i.e. positions 2_030 and 2_029.
    let record: VariantRecord = insertion(2_029, 2_030, Strand::Reverse, &TRANSCRIPT[10..30]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            2_030, Strand::Reverse, GraphOperationType::Upstream,
            2_049, Strand::Reverse, GraphOperationType::Downstream
        )
    ]);
}


/// An inverted piece on a reverse-strand transcript is traversed on the forward strand:
/// the mirror image of the forward-transcript inverted duplication.
#[test]
fn breakpoint_rescue_inverted_piece_on_reverse_strand_transcript_is_traversed_forward() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, (2_000..2_060).rev(), Strand::Reverse)
    ];
    let query: String = reverse_complement(&TRANSCRIPT[10..30]).to_string();
    let record: VariantRecord = insertion(2_029, 2_030, Strand::Reverse, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            2_030, Strand::Reverse, GraphOperationType::Upstream,
            2_030, Strand::Forward, GraphOperationType::Upstream
        ),
        breakpoint(
            2_029, Strand::Reverse, GraphOperationType::Downstream,
            2_049, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// Two exons, 1_000..1_030 and 5_000..5_030, so transcript offsets 29 and 30 are the exon
/// junction. An insertion at the junction copying the end of exon 1 is a tandem
/// duplication in transcript space even though its flanks are 4_000 bases apart, and a
/// copy of the start of exon 2 continues across the intron before jumping back.
#[test]
fn breakpoint_rescue_works_in_spliced_coordinates_across_an_intron() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, (1_000..1_030).chain(5_000..5_030), Strand::Forward)
    ];

    let record: VariantRecord = insertion(1_029, 5_000, Strand::Forward, &TRANSCRIPT[10..30]);
    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        )
    ]);

    let record: VariantRecord = insertion(1_029, 5_000, Strand::Forward, &TRANSCRIPT[30..50]);
    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            5_000, Strand::Forward, GraphOperationType::Upstream,
            5_019, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// The first transcript does not contain the flanks; the second does and is used.
#[test]
fn breakpoint_rescue_skips_transcripts_missing_a_flank() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 2_000..2_060, Strand::Forward),
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// Both transcripts hold the flanks, but the first has the duplicated bases masked out, so
/// nothing places on it and the rescue falls through to the second instead of giving up.
#[test]
fn breakpoint_rescue_falls_through_when_placement_fails_on_the_first_transcript() {
    let masked: String = format!("{}{}{}", &TRANSCRIPT[..10], "N".repeat(20), &TRANSCRIPT[30..]);
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(&masked, 1_000..1_060, Strand::Forward),
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_029, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// Flank lookup is per transcript. No transcripts, a missing left flank, or the two flanks
/// split across two transcripts all leave nothing to place on.
#[test]
fn breakpoint_rescue_needs_both_flanks_on_a_single_transcript() {
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);

    // No transcripts at all.
    assert!(retype(&record, &Vec::new(), 0.0, 3).is_none());

    // Left flank one base before the transcript start.
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_030..1_090, Strand::Forward)
    ];
    assert!(retype(&record, &transcripts, 0.0, 3).is_none());

    // Each transcript holds one flank: 1_029 ends the first, 1_030 starts the second.
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 970..1_030, Strand::Forward),
        build_transcript(TRANSCRIPT, 1_030..1_090, Strand::Forward)
    ];
    assert!(retype(&record, &transcripts, 0.0, 3).is_none());
}


/// Three exons; the read carries the middle one but the aligner reported it as an
/// insertion between the outer two. Every hop of the walk is continuation in transcript
/// space, so no breakpoint is produced and the rescue answers `None`, which leaves the
/// caller holding the original insertion record.
#[test]
fn breakpoint_rescue_insertion_that_is_the_skipped_exon_is_plain_continuation() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(
            TRANSCRIPT,
            (1_000..1_010).chain(3_000..3_020).chain(5_000..5_030),
            Strand::Forward
        )
    ];
    let record: VariantRecord = insertion(1_009, 5_000, Strand::Forward, &TRANSCRIPT[10..30]);

    assert!(retype(&record, &transcripts, 1.0, 3).is_none());
}


/// Eighteen of thirty bases placed is exactly 0.6: accepted at 0.6, rejected one hundredth
/// above it.
#[test]
fn breakpoint_rescue_placed_fraction_threshold_is_inclusive() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = format!("{}{}", &TRANSCRIPT[40..52], &TRANSCRIPT[2..20]);
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &query);

    assert!(retype(&record, &transcripts, 0.6, 1).is_some());
    assert!(retype(&record, &transcripts, 0.61, 1).is_none());
}


/// The flank index is keyed by position alone, so a transcript on another chromosome that
/// spans the same coordinates places the insertion and reports breakpoints on its own
/// chromosome. Ignored until the loop skips transcripts whose chromosome differs from the
/// insertion's.
#[test]
#[ignore]
fn breakpoint_rescue_ignores_transcripts_on_another_chromosome() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript_helper(1, TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Forward, &TRANSCRIPT[10..30]);

    let result: Option<Vec<VariantRecord>> = retype(&record, &transcripts, 0.9, 3);

    assert!(result.is_none(), "placed on the wrong chromosome: {:?}", result);
}

/// The tandem duplication of `breakpoint_rescue_tandem_duplication_yields_one_back_junction`,
/// carried by a read of the same molecule sequenced from its other end. The read is aligned on
/// the reverse strand and holds the copy reverse-complemented. Read in its own order it meets
/// base 30 first and runs back along the transcript, so it spells the one back junction with both
/// strands flipped, as a split alignment of it would, and not as an inverted copy.
#[test]
fn breakpoint_rescue_antisense_read_spells_the_tandem_duplication_with_both_strands_flipped() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let query: String = reverse_complement(&TRANSCRIPT[10..30]).to_string();
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Reverse, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Reverse, GraphOperationType::Upstream,
            1_029, Strand::Reverse, GraphOperationType::Downstream
        )
    ]);

    // The same on a reverse-strand transcript, read on the forward strand.
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, (2_000..2_060).rev(), Strand::Reverse)
    ];
    let record: VariantRecord = insertion(2_029, 2_030, Strand::Forward, &query);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            2_030, Strand::Forward, GraphOperationType::Upstream,
            2_049, Strand::Forward, GraphOperationType::Downstream
        )
    ]);
}


/// The inverted duplication of `breakpoint_rescue_inverted_duplication_flips_strand_across_the_piece`
/// carried by a read sequenced from the molecule's other end: it holds the transcript's own bases.
/// It is still an inversion (two junctions that change strand), met in the other order, and not
/// the tandem duplication its bases spell on the transcript's strand. The fold-back at base 29 has
/// both sides at one position, so its sides keep the read's order.
#[test]
fn breakpoint_rescue_antisense_read_keeps_an_inverted_duplication_inverted() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let record: VariantRecord = insertion(1_029, 1_030, Strand::Reverse, &TRANSCRIPT[10..30]);

    let records: Vec<VariantRecord> = retype(&record, &transcripts, 0.9, 3).unwrap();

    assert_eq!(records, vec![
        breakpoint(
            1_010, Strand::Forward, GraphOperationType::Upstream,
            1_030, Strand::Reverse, GraphOperationType::Upstream
        ),
        breakpoint(
            1_029, Strand::Forward, GraphOperationType::Downstream,
            1_029, Strand::Reverse, GraphOperationType::Downstream
        )
    ]);
}


/// A read aligned to bases 0..=29 that ends in a 20-base clip holding bases 5..25, as a second lap
/// gives. The clip is an insertion on the read's last bases. The read leaves the transcript at
/// base 29 and enters it at base 5, and ends inside the clip: it never comes back to base 30, so
/// the walk takes no closing hop. A read that starts with a clip holding bases 40..60 and is
/// aligned from base 30 on likewise only enters at base 30 from base 59.
#[test]
fn breakpoint_rescue_takes_no_hop_beyond_the_end_of_the_read_for_a_terminal_clip() {
    let transcripts: Vec<ReferenceTranscriptSequence> = vec![
        build_transcript(TRANSCRIPT, 1_000..1_060, Strand::Forward)
    ];
    let clip = |read_position_1: u32, read_position_2: u32, sequence: &str| -> VariantRecord {
        VariantRecord::new(
            READ_ID,
            read_position_1,
            read_position_2,
            GraphOperation::new(
                CHROMOSOME, 1_029, Strand::Forward, GraphOperationType::Downstream,
                CHROMOSOME, 1_030, Strand::Forward, GraphOperationType::Upstream,
                sequence.into(),
                VariantType::Insertion
            )
        )
    };
    let operations = |records: Vec<VariantRecord>| -> Vec<Box<str>> {
        records.iter().map(|record| record.get_graph_operation().as_boxed_str()).collect()
    };

    // Trailing clip: the read's last 20 bases of 1,000.
    let records: Vec<VariantRecord> = retype(&clip(980, 999, &TRANSCRIPT[5..25]), &transcripts, 0.9, 3).unwrap();
    assert_eq!(operations(records), vec![
        breakpoint(1_005, Strand::Forward, GraphOperationType::Upstream, 1_029, Strand::Forward, GraphOperationType::Downstream)
            .get_graph_operation()
            .as_boxed_str()
    ]);

    // Leading clip: the read's first 20 bases.
    let records: Vec<VariantRecord> = retype(&clip(0, 19, &TRANSCRIPT[40..60]), &transcripts, 0.9, 3).unwrap();
    assert_eq!(operations(records), vec![
        breakpoint(1_030, Strand::Forward, GraphOperationType::Upstream, 1_059, Strand::Forward, GraphOperationType::Downstream)
            .get_graph_operation()
            .as_boxed_str()
    ]);

    // The same bases inside the read take both hops.
    let records: Vec<VariantRecord> = retype(&clip(100, 119, &TRANSCRIPT[5..25]), &transcripts, 0.9, 3).unwrap();
    assert_eq!(records.len(), 2);
}
