use super::*;

fn transcript(id: &str, reverse: bool, exons: &[(u32, u32)]) -> Transcript {
    let strand = if reverse { Strand::Reverse } else { Strand::Forward };
    let mut tx = Transcript::new(
        "GENE", id, "test", "chr1", exons[0].0, exons.last().unwrap().1,
        strand.clone(), 1, id, "protein_coding", "1", HashSet::new()
    );
    let mut ordered = exons.to_vec();
    if reverse { ordered.reverse(); }
    for (index, &(start, end)) in ordered.iter().enumerate() {
        tx.add_exon(Exon::new(
            "GENE", id, &format!("EX{index}"), "test", "chr1", start, end,
            strand.clone(), 1, index as u16 + 1
        ));
    }
    tx
}

#[test]
fn trims_internal_exon_overhangs_at_both_ends_on_both_strands() {
    for reverse in [false, true] {
        for (sequence, cigar, cs, exons, expected) in [
            ("AAAAACCCCCCCC", "5M10N8M", ":5~gt10ag:8",
                [(81, 105), (116, 120), (131, 140)], "AAAAACCCCC"),
            ("AAAAAGGGGGCCCCC", "10M10N5M", ":10~gt10ag:5",
                [(81, 90), (106, 110), (121, 140)], "GGGGGCCCCC")
        ] {
            let read = model(sequence, cigar, cs, reverse);
            let annotations = TranscriptAnnotations(vec![transcript("tx", reverse, &exons)]);
            let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
            assert_eq!(corrected.base_quality_scores, vec![40; expected.len()]);
        }
    }
}

#[test]
fn true_transcript_ends_and_unmatched_introns_keep_the_overhang() {
    for reverse in [false, true] {
        let read = model("AAAAACCCCCCCC", "5M10N8M", ":5~gt10ag:8", reverse);
        for annotations in [
            TranscriptAnnotations(vec![transcript("terminal", reverse, &[(81, 105), (116, 120)])]),
            TranscriptAnnotations(vec![transcript("wrong_intron", reverse, &[(81, 104), (116, 120), (131, 140)])]),
            TranscriptAnnotations(vec![transcript("wrong_strand", !reverse, &[(81, 105), (116, 120), (131, 140)])]),
            TranscriptAnnotations::default()
        ] {
            let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
            assert_eq!(&*corrected.sequence, read.get_read_sequence());
        }
    }
}

#[test]
fn alternative_isoforms_choose_the_outermost_boundary_or_veto_trimming() {
    for reverse in [false, true] {
        let read = model("AAAAACCCCCCCC", "5M10N8M", ":5~gt10ag:8", reverse);
        for (alternative, expected) in [
            (vec![(81, 105), (116, 122), (131, 140)], "AAAAACCCCCCC"),
            (vec![(81, 105), (116, 124), (131, 140)], "AAAAACCCCCCCC"),
            (vec![(81, 105), (116, 120)], "AAAAACCCCCCCC")
        ] {
            let annotations = TranscriptAnnotations(vec![
                transcript("short", reverse, &[(81, 105), (116, 120), (131, 140)]),
                transcript("alternative", reverse, &alternative)
            ]);
            let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
        }
    }
}

// A trusted variant in the span a transcript end would remove vetoes the trim, whether it
// names a position inside the span or only covers it, and whether it is small or a structural
// anchor. The trusted allele is then applied, so the veto shows as the overhang surviving: an
// SNV or MNV the read already carries leaves it as sequenced, an insertion lands inside it, and
// a deletion removes only its own span from it. Without the veto the read ends at 120.
#[test]
fn confident_variants_protect_the_removed_span_including_structural_anchors() {
    for reverse in [false, true] {
        let read = model("AAAAACCCCCCCC", "5M10N8M", ":5~gt10ag:8", reverse);
        let annotations = TranscriptAnnotations(vec![transcript("tx", reverse, &[(81, 105), (116, 120), (131, 140)])]);
        for (protected, expected) in [
            (operation(121, 123, "C", VariantType::SingleNucleotideVariant), "AAAAACCCCCCCC"),
            // Neither MNV endpoint lies in 121..123, but its reference span does.
            (operation(120, 124, "CCC", VariantType::MultiNucleotideVariant), "AAAAACCCCCCCC"),
            (operation(120, 121, "G", VariantType::Insertion), "AAAAACCCCCGCCC"),
            (operation(122, 124, "", VariantType::Deletion), "AAAAACCCCCCC"),
            (GraphOperation::new(0, 122, Strand::Forward, GraphOperationType::Downstream,
                1, 200, Strand::Forward, GraphOperationType::Upstream, "".into(), VariantType::FusionGene),
                "AAAAACCCCCCCC")
        ] {
            let corrected = correct_annotated(&read, &HashSet::from([protected]), &annotations, 50, 7);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
        }
        let other_chromosome = GraphOperation::new(
            1, 121, Strand::Forward, GraphOperationType::Downstream,
            1, 123, Strand::Forward, GraphOperationType::Upstream, "C".into(), VariantType::SingleNucleotideVariant
        );
        let corrected = correct_annotated(&read, &HashSet::from([other_chromosome]), &annotations, 50, 7);
        assert_eq!(&*corrected.sequence, oriented("AAAAACCCCC", reverse));
    }
}

#[test]
fn transcript_trims_take_their_own_clip_and_keep_the_other() {
    for reverse in [false, true] {
        let read = model("GGAAAAACCCCCCCCTT", "2S5M10N8M2S", ":5~gt10ag:8", reverse);
        let annotations = TranscriptAnnotations(vec![transcript("tx", reverse, &[(81, 105), (116, 120), (131, 140)])]);
        let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
        assert_eq!(&*corrected.sequence, oriented("GGAAAAACCCCC", reverse));

        let protected_clip = HashSet::from([operation(123, 124, "TT", VariantType::Insertion)]);
        let corrected = correct_annotated(&read, &protected_clip, &annotations, 50, 7);
        assert_eq!(&*corrected.sequence, oriented("GGAAAAACCCCCCCCTT", reverse));
    }
}

#[test]
fn deletions_crossing_a_trim_boundary_are_not_restored() {
    for reverse in [false, true] {
        let read = model("AAAAACCCCGG", "5M10N4M2D2M", ":5~gt10ag:4-aa:2", reverse);
        let annotations = TranscriptAnnotations(vec![transcript("tx", reverse, &[(81, 105), (116, 120), (131, 140)])]);
        let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
        assert_eq!(&*corrected.sequence, oriented("AAAAACCCC", reverse));
        assert_eq!(corrected.base_quality_scores, [40; 9]);
    }
}

#[test]
fn corrections_inside_the_kept_span_still_receive_the_requested_quality() {
    for reverse in [false, true] {
        let read = model("AAAAACCTCCCCC", "5M10N8M", ":5~gt10ag:2*ct:5", reverse);
        let annotations = TranscriptAnnotations(vec![transcript("tx", reverse, &[(81, 105), (116, 120), (131, 140)])]);
        let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
        assert_eq!(&*corrected.sequence, oriented("AAAAACCCCC", reverse));
        let mut quality = vec![40, 40, 40, 40, 40, 40, 40, 7, 40, 40];
        if reverse { quality.reverse(); }
        assert_eq!(corrected.base_quality_scores, quality);
    }
}

#[test]
fn single_block_reads_are_not_trimmed_to_annotation() {
    let read = model("AAAAACCCCCCCC", "13M", ":13", false);
    let annotations = TranscriptAnnotations(vec![transcript("tx", false, &[(81, 105), (116, 120), (131, 140)])]);
    let corrected = correct_annotated(&read, &HashSet::new(), &annotations, 50, 7);
    assert_eq!(&*corrected.sequence, read.get_read_sequence());
}
