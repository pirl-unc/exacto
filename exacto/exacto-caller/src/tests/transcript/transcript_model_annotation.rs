use super::*;


/// Runs of read positions appended in read order. A run that starts at the position after the
/// last one and holds the same annotation extends it. A run after a gap, or with another
/// annotation, is a run of its own.
#[test]
fn push_bases_returns_matches() {
    let mut annotation: TranscriptModelAnnotation = TranscriptModelAnnotation::new();
    annotation.push_bases(0, 9, TranscriptModelBaseAnnotation {
        context: AlignmentModelBaseContext::Exonic,
        reference_gene_id: Some("gene".into()),
        reference_transcript_id: Some("transcript".into()),
        reference_exon_id: Some("exon_1".into())
    });
    annotation.push_bases(10, 19, TranscriptModelBaseAnnotation {
        context: AlignmentModelBaseContext::Exonic,
        reference_gene_id: Some("gene".into()),
        reference_transcript_id: Some("transcript".into()),
        reference_exon_id: Some("exon_1".into())
    });
    annotation.push_bases(25, 29, TranscriptModelBaseAnnotation {
        context: AlignmentModelBaseContext::Exonic,
        reference_gene_id: Some("gene".into()),
        reference_transcript_id: Some("transcript".into()),
        reference_exon_id: Some("exon_1".into())
    });
    annotation.push_bases(30, 39, TranscriptModelBaseAnnotation {
        context: AlignmentModelBaseContext::Exonic,
        reference_gene_id: Some("gene".into()),
        reference_transcript_id: Some("transcript".into()),
        reference_exon_id: Some("exon_2".into())
    });
    annotation.push_bases(40, 40, TranscriptModelBaseAnnotation {
        context: AlignmentModelBaseContext::Intergenic,
        reference_gene_id: None,
        reference_transcript_id: None,
        reference_exon_id: None
    });

    // (first read position, last read position, context, exon ID)
    let runs: Vec<(u32, u32, AlignmentModelBaseContext, Option<&str>)> = annotation
        .bases
        .iter()
        .map(|run| (run.read_start, run.read_end, run.annotation.context.clone(), run.annotation.reference_exon_id.as_deref()))
        .collect();
    assert_eq!(
        runs,
        vec![
            (0, 19, AlignmentModelBaseContext::Exonic, Some("exon_1")),
            (25, 29, AlignmentModelBaseContext::Exonic, Some("exon_1")),
            (30, 39, AlignmentModelBaseContext::Exonic, Some("exon_2")),
            (40, 40, AlignmentModelBaseContext::Intergenic, None)
        ]
    );
    assert_eq!(
        annotation.get_read_positions().collect::<Vec<u32>>(),
        (0..20).chain(25..41).collect::<Vec<u32>>()
    );
    assert_eq!(annotation.get_base(19).unwrap().reference_exon_id.as_deref(), Some("exon_1"));
    assert!(annotation.get_base(20).is_none());
    assert_eq!(annotation.get_base(30).unwrap().reference_exon_id.as_deref(), Some("exon_2"));

    // The same annotation, one read position at a time.
    let mut annotation_per_base: TranscriptModelAnnotation = TranscriptModelAnnotation::new();
    for read_position in annotation.get_read_positions().collect::<Vec<u32>>() {
        annotation_per_base.set_base(read_position, annotation.get_base(read_position).unwrap().clone());
    }
    let runs_per_base: Vec<(u32, u32, AlignmentModelBaseContext, Option<&str>)> = annotation_per_base
        .bases
        .iter()
        .map(|run| (run.read_start, run.read_end, run.annotation.context.clone(), run.annotation.reference_exon_id.as_deref()))
        .collect();
    assert_eq!(runs_per_base, runs);
}


/// Runs of skipped reference positions added to an event in any order. Runs that follow one
/// another with the same annotation become one run, with the nucleotides in position order.
/// A run with another annotation stays a run of its own. The event is the same whichever of
/// its two read positions is named first.
#[test]
fn add_skipped_reference_run_returns_matches() {
    let mut annotation: TranscriptModelAnnotation = TranscriptModelAnnotation::new();
    // (read position 1, read position 2, first position, last position, exon ID, nucleotides)
    for (read_position_1, read_position_2, reference_start, reference_end, reference_exon_id, sequence) in [
        (10u32, 11u32, 100u32, 104u32, "exon_1", "ACGTA"),
        (10, 11, 110, 112, "exon_1", "GGG"),
        (11, 10, 105, 109, "exon_1", "CCCCC"),
        (10, 11, 113, 115, "exon_2", "TTT")
    ] {
        annotation.add_skipped_reference_run(read_position_1, read_position_2, SkippedReferenceRun {
            reference_chromosome_id: 0,
            reference_start,
            reference_end,
            reference_strand: Strand::Forward,
            reference_gene_id: Some("gene".into()),
            reference_transcript_id: Some("transcript".into()),
            reference_exon_id: Some(reference_exon_id.into()),
            sequence: sequence.into()
        });
    }

    // (first position, last position, exon ID, nucleotides)
    let runs: Vec<(u32, u32, Option<&str>, &str)> = annotation
        .get_event(10, 11)
        .unwrap()
        .get_skipped_runs()
        .iter()
        .map(|run| (run.reference_start, run.reference_end, run.reference_exon_id.as_deref(), &*run.sequence))
        .collect();
    assert_eq!(
        runs,
        vec![
            (100, 112, Some("exon_1"), "ACGTACCCCCGGG"),
            (113, 115, Some("exon_2"), "TTT")
        ]
    );
    assert_eq!(annotation.events.len(), 1);

    // The same positions, one reference base at a time and out of order.
    let mut annotation_per_base: TranscriptModelAnnotation = TranscriptModelAnnotation::new();
    for reference_base in annotation
        .get_event(10, 11)
        .unwrap()
        .get_skipped_reference_bases()
        .into_iter()
        .flatten()
        .rev()
    {
        annotation_per_base.add_skipped_reference_base(10, 11, reference_base);
    }
    let runs_per_base: Vec<(u32, u32, Option<&str>, &str)> = annotation_per_base
        .get_event(10, 11)
        .unwrap()
        .get_skipped_runs()
        .iter()
        .map(|run| (run.reference_start, run.reference_end, run.reference_exon_id.as_deref(), &*run.sequence))
        .collect();
    assert_eq!(runs_per_base, runs);
}
