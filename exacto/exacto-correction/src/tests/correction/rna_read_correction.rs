use std::sync::Arc;
use exacto_caller::alignment::alignment_model_editor::AlignmentModelEditor;
use exacto_core::prelude::{get_bam_fastx_base_quality_scores, get_bam_fastx_read_sequence, reverse_complement};
use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write;
use super::*;

/// Build real BAM-backed models from small, explicit CIGAR/cs examples. SAM sequence
/// is in reference orientation; AlignmentModel receives the original read orientation.
fn model(sequence: &str, cigar: &str, cs: &str, reverse: bool) -> AlignmentModel {
    let sam_text = format!(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-7\t{}\tchr1\t101\t60\t{}\t*\t0\t0\t{}\t{}\tcs:Z:{}\n",
        if reverse { 16 } else { 0 }, cigar, sequence, "I".repeat(sequence.len()), cs
    );
    let mut reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = reader.read_header().unwrap();
    let record = reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records = vec![record];
    let sequence = get_bam_fastx_read_sequence(&records);
    let quality = get_bam_fastx_base_quality_scores(&records);
    AlignmentModel::new(7, &sequence, &quality, &records.into_iter().map(Arc::new).collect::<Vec<_>>())
}

fn operation(left: u32, right: u32, sequence: &str, kind: VariantType) -> GraphOperation {
    GraphOperation::new(0, left, Strand::Forward, GraphOperationType::Downstream,
        0, right, Strand::Forward, GraphOperationType::Upstream, sequence.into(), kind)
}

// Check the public output and the complete input model, including base metadata,
// event keys and the endpoint index, after every fixture correction.
fn correct(read: &AlignmentModel, confident: &HashSet<GraphOperation>, max_len: usize, quality: u8) -> CorrectedRNARead {
    correct_annotated(read, confident, &TranscriptAnnotations::default(), max_len, quality)
}

fn correct_annotated(read: &AlignmentModel, confident: &HashSet<GraphOperation>, annotator: &TranscriptAnnotations, max_len: usize, quality: u8) -> CorrectedRNARead {
    let before = format!("{read:?}");
    let chromosome_names: BiMap<Box<str>, u16> = BiMap::from_iter([("chr1".into(), 0)]);
    let confident: HashSet<&GraphOperation> = confident.iter().collect();
    let corrected = correct_rna_read(42, "read-7", read, &confident, &chromosome_names, Some(annotator), max_len, quality);
    assert_eq!(format!("{read:?}"), before, "correction must not change the model");
    assert_eq!(corrected.sequence.len(), corrected.base_quality_scores.len());
    assert_eq!(corrected.cluster_id, 42);
    assert_eq!(&*corrected.read_name, "read-7");
    corrected
}

#[derive(Default)]
struct TranscriptAnnotations(Vec<Transcript>);

impl GeneAnnotator for TranscriptAnnotations {
    fn get_assembly(&self) -> &str { "test" }
    fn get_version(&self) -> &str { "test" }
    fn get_gene_ids_at_locus(&self, _: &str, _: u32) -> Vec<Box<str>> { Vec::new() }
    fn get_gene_ids_overlapping_region(&self, _: &str, _: u32, _: u32) -> Vec<Box<str>> { Vec::new() }
    fn get_transcript_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>> {
        self.0.iter().filter(|tx| &*tx.chromosome == chromosome && tx.start <= end && start <= tx.end)
            .map(|tx| tx.transcript_id.clone()).collect()
    }
    fn get_exon_ids_overlapping_region(&self, _: &str, _: u32, _: u32) -> Vec<Box<str>> { Vec::new() }
    fn get_gene(&self, _: &str) -> Option<&Gene> { None }
    fn get_genes(&self) -> Vec<&Gene> { Vec::new() }
    fn get_transcript(&self, id: &str) -> Option<&Transcript> {
        self.0.iter().find(|tx| &*tx.transcript_id == id)
    }
    fn get_transcripts(&self) -> Vec<&Transcript> { self.0.iter().collect() }
    fn get_exon(&self, tx: &str, exon: &str) -> Option<&Exon> { self.get_transcript(tx)?.get_exon(exon) }
    fn get_exons(&self) -> Vec<&Exon> { self.0.iter().flat_map(Transcript::get_sorted_exons).collect() }
    fn rank_transcripts<'a>(&self, transcripts: Vec<&'a Transcript>) -> Vec<&'a Transcript> { transcripts }
}

#[path = "rna_transcript_end_trimming.rs"]
mod transcript_ends;

fn oriented(sequence: &str, reverse: bool) -> String {
    if reverse { reverse_complement(sequence).to_string() } else { sequence.to_string() }
}

#[test]
fn unsupported_snv_insertion_and_deletion_are_corrected_on_both_strands() {
    for reverse in [false, true] {
        let read = model("ACCGTG", "3M1I1M2D1M", ":1*ac:1+g:1-aa:1", reverse);
        let corrected = correct(&read, &HashSet::new(), 50, 0);
        assert_eq!(&*corrected.sequence, oriented("AACTAAG", reverse));
        assert_eq!(corrected.base_quality_scores.iter().filter(|&&q| q == 0).count(), 3);
    }
}

#[test]
fn confident_alleles_are_preserved_even_when_the_read_is_on_the_other_strand() {
    let alleles = [
        operation(101, 103, "C", VariantType::SingleNucleotideVariant),
        operation(103, 104, "G", VariantType::Insertion),
        operation(104, 107, "", VariantType::Deletion)
    ];
    for reverse in [false, true] {
        let read = model("ACCGTG", "3M1I1M2D1M", ":1*ac:1+g:1-aa:1", reverse);
        for (allele, expected) in alleles.iter().zip(["ACCTAAG", "AACGTAAG", "AACTG"]) {
            let corrected = correct(&read, &HashSet::from([allele.clone()]), 50, 0);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
        }
        let corrected = correct(&read, &HashSet::from(alleles.clone()), 50, 0);
        assert_eq!(&*corrected.sequence, read.get_read_sequence());
        assert_eq!(corrected.base_quality_scores, [40; 6]);
    }
}

#[test]
fn neighboring_error_does_not_erase_a_confident_snv() {
    for reverse in [false, true] {
        let read = model("ACTA", "4M", ":1*ac*gt:1", reverse);
        let confident = HashSet::from([operation(101, 103, "C", VariantType::SingleNucleotideVariant)]);
        let corrected = correct(&read, &confident, 50, 0);
        assert_eq!(&*corrected.sequence, oriented("ACGA", reverse));
    }
}

// The cluster's spelling replaces the read's own at every trusted site, even where the read
// showed something else: the read's C becomes the called T, its inserted G becomes the called A,
// and of its two-base deletion only the called one-base span stays deleted. The corrected read
// is "reference plus the cluster's alleles", and the replaced bases carry the requested quality.
#[test]
fn trusted_alleles_replace_the_reads_own_spelling_at_each_site() {
    let confident = HashSet::from([
        operation(101, 103, "T", VariantType::SingleNucleotideVariant),
        operation(103, 104, "A", VariantType::Insertion),
        operation(104, 106, "", VariantType::Deletion)
    ]);
    let read = model("ACCGTG", "3M1I1M2D1M", ":1*ac:1+g:1-aa:1", false);
    let corrected = correct(&read, &confident, 50, 0);
    assert_eq!(&*corrected.sequence, "ATCATAG");
    assert_eq!(corrected.base_quality_scores, [40, 0, 40, 0, 40, 0, 40]);
}

#[test]
fn deletion_payload_is_emitted_even_when_its_boundary_base_is_skipped() {
    // Test both CIGAR orders and both strands. In at least one read orientation the
    // deletion boundary falls on an insertion base that will not be emitted.
    for (cigar, cs) in [("2M1I2D2M", ":2+g-aa:2"), ("2M2D1I2M", ":2-aa+g:2")] {
        for reverse in [false, true] {
            let read = model("ACGTT", cigar, cs, reverse);
            let corrected = correct(&read, &HashSet::new(), 50, 7);
            assert_eq!(&*corrected.sequence, oriented("ACAATT", reverse));
            let mut expected_quality = vec![40, 40, 7, 7, 40, 40];
            if reverse { expected_quality.reverse(); }
            assert_eq!(corrected.base_quality_scores, expected_quality);
        }
    }
}

#[test]
fn event_length_limit_is_respected_including_zero() {
    let read = model("ACCGTG", "3M1I1M2D1M", ":1*ac:1+g:1-aa:1", false);
    let corrected = correct(&read, &HashSet::new(), 1, 0);
    assert_eq!(&*corrected.sequence, "AACTG");
    let corrected = correct(&read, &HashSet::new(), 0, 0);
    assert_eq!(&*corrected.sequence, "ACCGTG");
    assert_eq!(corrected.base_quality_scores, [40; 6]);
}

// An uncalled terminal clip is kept as sequenced. A clip the cluster called as an insertion at
// its anchor is kept with the read's own qualities when the spelling matches, and replaced by the
// cluster's spelling at the requested quality when it does not.
#[test]
fn terminal_clips_are_kept_or_replaced_by_the_trusted_insertion_on_both_strands() {
    for reverse in [false, true] {
        let read = model("GGACCTT", "2S3M2S", ":3", reverse);
        for (confident, expected, expected_quality) in [
            (HashSet::new(), "GGACCTT", vec![40; 7]),
            (HashSet::from([operation(100, 101, "GG", VariantType::Insertion)]), "GGACCTT", vec![40; 7]),
            (HashSet::from([operation(103, 104, "TT", VariantType::Insertion)]), "GGACCTT", vec![40; 7]),
            (HashSet::from([
                operation(100, 101, "GG", VariantType::Insertion),
                operation(103, 104, "TT", VariantType::Insertion)
            ]), "GGACCTT", vec![40; 7]),
            (HashSet::from([operation(100, 101, "AA", VariantType::Insertion)]), "AAACCTT", vec![0, 0, 40, 40, 40, 40, 40])
        ] {
            let corrected = correct(&read, &confident, 50, 0);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
            let mut expected_quality = expected_quality;
            if reverse { expected_quality.reverse(); }
            assert_eq!(corrected.base_quality_scores, expected_quality);
        }
    }
}

// The event length cap governs internal insertions, not terminal clips: an uncalled clip beyond
// the cap is kept like any other, and a trusted insertion beyond the cap still replaces its
// clip, base for base, at the requested quality. An uncalled fusion arm is kept as sequenced
// and a called one is written in the cluster's spelling, at any scale.
#[test]
fn terminal_clips_and_their_trusted_replacements_ignore_the_event_length_cap() {
    let clip = "G".repeat(60);
    let called = "T".repeat(60);
    for reverse in [false, true] {
        let read = model(&format!("ACC{clip}"), "3M60S", ":3", reverse);
        let corrected = correct(&read, &HashSet::new(), 50, 7);
        assert_eq!(&*corrected.sequence, read.get_read_sequence());
        assert_eq!(corrected.base_quality_scores, [40; 63]);

        let confident = HashSet::from([operation(103, 104, &called, VariantType::Insertion)]);
        let corrected = correct(&read, &confident, 50, 7);
        assert_eq!(&*corrected.sequence, oriented(&format!("ACC{called}"), reverse));
        let mut expected_quality = [vec![40; 3], vec![7; 60]].concat();
        if reverse { expected_quality.reverse(); }
        assert_eq!(corrected.base_quality_scores, expected_quality);
    }
}

#[test]
fn empty_and_unplaced_reads_keep_the_original_sequence_and_quality() {
    for sequence in ["", "ACGT"] {
        let read = AlignmentModel::from_read(7, sequence, &vec![31; sequence.len()]);
        let corrected = correct(&read, &HashSet::new(), 50, 0);
        assert_eq!(&*corrected.sequence, sequence);
        assert_eq!(corrected.base_quality_scores, vec![31; sequence.len()]);
    }
}

#[test]
fn clipping_and_aligned_corrections_share_original_coordinates() {
    for reverse in [false, true] {
        let read = model("GGACCGTGTT", "2S3M1I1M2D1M2S", ":1*ac:1+g:1-aa:1", reverse);
        let corrected = correct(&read, &HashSet::new(), 50, 13);
        assert_eq!(&*corrected.sequence, oriented("GGAACTAAGTT", reverse));
        let mut expected_quality = vec![40, 40, 40, 13, 40, 40, 13, 13, 40, 40, 40];
        if reverse { expected_quality.reverse(); }
        assert_eq!(corrected.base_quality_scores, expected_quality);
    }
}

#[test]
fn splice_events_and_their_inserted_sequence_are_preserved() {
    for (sequence, cigar, cs) in [
        ("ACGT", "2M10N2M", ":2~gt10ag:2"),
        ("ACGTT", "2M1I10N2M", ":2+g~gt10ag:2")
    ] {
        for reverse in [false, true] {
            let read = model(sequence, cigar, cs, reverse);
            let corrected = correct(&read, &HashSet::new(), 50, 0);
            assert_eq!(&*corrected.sequence, read.get_read_sequence());
            assert!(read.is_spliced());
        }
    }
}

// A deletion on each side of an intron, `-ac~gt30ag-tt`, must never fill the intron. The
// constructor files one Splicing event for those flanks (a splice supersedes a deletion at the
// same key), so no deletion record exists and neither the intron nor the two deletions are
// restored. The 30 bp intron keeps the whole span under the length cap, so nothing but the
// rule itself protects the junction.
#[test]
fn deletions_flanking_a_splice_never_fill_the_intron() {
    for reverse in [false, true] {
        let read = model("ACGT", "2=2D30N2D2=", ":2-ac~gt30ag-tt:2", reverse);
        assert!(read.is_spliced());
        assert_eq!(read.num_events(), 1, "the splice must be the only event on those flanks");
        let corrected = correct(&read, &HashSet::new(), 50, 7);
        assert_eq!(&*corrected.sequence, oriented("ACGT", reverse), "bases were invented across the junction");
        assert_eq!(corrected.base_quality_scores, [40; 4]);
    }
}

#[test]
fn confident_mnv_in_reverse_orientation_preserves_forward_and_reverse_reads() {
    let confident = HashSet::from([GraphOperation::new(
        0, 101, Strand::Reverse, GraphOperationType::Downstream,
        0, 104, Strand::Reverse, GraphOperationType::Upstream,
        "AG".into(), VariantType::MultiNucleotideVariant
    )]);
    for reverse in [false, true] {
        let read = model("ACTA", "4M", ":1*ac*gt:1", reverse);
        let corrected = correct(&read, &confident, 50, 0);
        assert_eq!(&*corrected.sequence, read.get_read_sequence());
        assert_eq!(corrected.base_quality_scores, [40; 4]);
    }
}

#[test]
fn shared_deletion_endpoints_and_terminal_mismatches_are_corrected() {
    for reverse in [false, true] {
        let read = model("CGT", "1M1D1M1D1M", "*ac-t:1-a*gt", reverse);
        assert_eq!(read.num_events(), 2);
        let corrected = correct(&read, &HashSet::new(), 50, 0);
        assert_eq!(&*corrected.sequence, oriented("ATGAG", reverse));
        assert_eq!(corrected.base_quality_scores, [0, 0, 40, 0, 0]);
    }
}

#[test]
fn entirely_soft_clipped_and_internal_unplaced_bases_are_retained() {
    let mut read = AlignmentModel::from_read(7, "", &[]);
    AlignmentModelEditor::new(&mut read).insert_unplaced(0, &[Nucleotide::A, Nucleotide::C], 19);
    let corrected = correct(&read, &HashSet::new(), 50, 0);
    assert_eq!(&*corrected.sequence, "AC");
    assert_eq!(corrected.base_quality_scores, [19, 19]);

    let mut read = model("ACGT", "4M", ":4", false);
    AlignmentModelEditor::new(&mut read).insert_unplaced(2, &[Nucleotide::T], 19);
    let corrected = correct(&read, &HashSet::new(), 50, 0);
    assert_eq!(&*corrected.sequence, "ACTGT");
    assert_eq!(corrected.base_quality_scores, [40, 40, 19, 40, 40]);
}

#[test]
fn orphan_deletion_payload_without_an_event_is_not_restored() {
    let mut read = model("ACGT", "2M2D2M", ":2-aa:2", false);
    AlignmentModelEditor::new(&mut read).remove(2..3);
    assert!(read.get_events().is_empty());
    assert_eq!(read.get_base(1).get_deletion_read_position(), Some(2));
    let corrected = correct(&read, &HashSet::new(), 50, 0);
    assert_eq!(&*corrected.sequence, "ACT");
}

#[test]
fn restored_deletion_stays_on_the_correct_side_of_a_confident_insertion() {
    for (cigar, cs, anchor, expected) in [
        ("2M1I2D2M", ":2+g-aa:2", 102, "ACGAATT"),
        ("2M2D1I2M", ":2-aa+g:2", 104, "ACAAGTT")
    ] {
        for reverse in [false, true] {
            let read = model("ACGTT", cigar, cs, reverse);
            let confident = HashSet::from([operation(anchor, anchor + 1, "G", VariantType::Insertion)]);
            let corrected = correct(&read, &confident, 50, 7);
            assert_eq!(&*corrected.sequence, oriented(expected, reverse));
        }
    }
}

// A called insertion or deletion, and reads whose aligner placed the same event at another copy
// of its tandem repeat, before or after the call. Below the length cap the read's own copy is an
// ordinary error and the call is written at its anchor. Above it the read keeps its own copy and
// is not given the call as well. Either way the event is in the corrected read once. The repeat
// unit is rna-004's; the reference is P + U x3 + S for the insertion, P + U x(3 + copies) + S for
// the deletion, from position 101, so the copies start after 108 or after 144.
#[test]
fn a_called_event_placed_elsewhere_in_its_repeat_is_applied_once() {
    let (p, s, unit) = ("GATTACAG", "TGACCTAG", "CCCATCCGCCTG");
    for reverse in [false, true] {
        for copies in [4usize, 5] {
            let event = unit.repeat(copies);
            let length = event.len() as u32;

            let read = format!("{p}{}{s}", unit.repeat(3 + copies));
            for called in [108, 144] {
                let confident = HashSet::from([operation(called, called + 1, &event, VariantType::Insertion)]);
                for (cigar, cs) in [
                    (format!("8M{length}I44M"), format!(":8+{}:44", event.to_lowercase())),
                    (format!("44M{length}I8M"), format!(":44+{}:8", event.to_lowercase()))
                ] {
                    let corrected = correct(&model(&read, &cigar, &cs, reverse), &confident, 50, 7);
                    assert_eq!(&*corrected.sequence, oriented(&read, reverse), "insertion of {length} called after {called}, {cigar}");
                }
            }

            let read = format!("{p}{}{s}", unit.repeat(3));
            for called in [108, 144] {
                let confident = HashSet::from([operation(called, called + length + 1, "", VariantType::Deletion)]);
                for (cigar, cs) in [
                    (format!("8M{length}D44M"), format!(":8-{}:44", event.to_lowercase())),
                    (format!("44M{length}D8M"), format!(":44-{}:8", event.to_lowercase()))
                ] {
                    let corrected = correct(&model(&read, &cigar, &cs, reverse), &confident, 50, 7);
                    assert_eq!(&*corrected.sequence, oriented(&read, reverse), "deletion of {length} called after {called}, {cigar}");
                }
            }
        }
    }
}

// A long copy at the call's own anchor is still written in the cluster's spelling: here the read
// lost one base of the 60-base insertion.
#[test]
fn a_long_copy_at_the_called_anchor_takes_the_called_spelling() {
    let (p, s, unit) = ("GATTACAG", "TGACCTAG", "CCCATCCGCCTG");
    let called = unit.repeat(5);
    let confident = HashSet::from([operation(108, 109, &called, VariantType::Insertion)]);
    for reverse in [false, true] {
        let read = format!("{p}{}{}{s}", &called[1..], unit.repeat(3));
        let cs = format!(":8+{}:44", called[1..].to_lowercase());
        let corrected = correct(&model(&read, "8M59I44M", &cs, reverse), &confident, 50, 7);
        assert_eq!(&*corrected.sequence, oriented(&format!("{p}{called}{}{s}", unit.repeat(3)), reverse));
    }
}

// A read that ends inside a called insertion holds its copy as a terminal clip: the aligner
// places two bases of the repeat and clips the rest. The read keeps its clip and is not given the
// call as well (rna-004 read 183).
#[test]
fn a_read_ending_inside_a_called_insertion_keeps_its_clip_and_is_not_given_the_call() {
    let (p, unit) = ("GATTACAG", "CCCATCCGCCTG");
    let confident = HashSet::from([operation(108, 109, &unit.repeat(5), VariantType::Insertion)]);
    for reverse in [false, true] {
        let read = format!("{p}{}", unit.repeat(6));
        let corrected = correct(&model(&read, "10M70S", ":10", reverse), &confident, 50, 7);
        assert_eq!(&*corrected.sequence, oriented(&read, reverse));
    }
}

// A long insertion farther from a called insertion than the call's length is not the call's
// copy: it is kept above the cap, and the call is still written at its own anchor.
#[test]
fn a_long_insertion_away_from_the_called_insertion_is_kept() {
    let (p, s, unit) = ("GATTACAG", "TGACCTAG", "CCCATCCGCCTG");
    let called = unit.repeat(5);
    let other = "ACGTTGCA".repeat(7) + "ACGT";
    let spacer = "ACGT".repeat(15);
    let confident = HashSet::from([operation(108, 109, &called, VariantType::Insertion)]);
    for reverse in [false, true] {
        // The read's insertion follows position 212, 104 bases from the call's anchor.
        let read = format!("{p}{}{s}{spacer}{other}TTGCAA", unit.repeat(3));
        let cs = format!(":112+{}:6", other.to_lowercase());
        let corrected = correct(&model(&read, "112M60I6M", &cs, reverse), &confident, 50, 7);
        let expected = format!("{p}{called}{}{s}{spacer}{other}TTGCAA", unit.repeat(3));
        assert_eq!(&*corrected.sequence, oriented(&expected, reverse));
    }
}
