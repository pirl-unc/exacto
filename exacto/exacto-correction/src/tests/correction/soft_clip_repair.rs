use exacto_core::prelude::*;
use std::collections::HashSet;
use tempfile::{tempdir, TempDir};

use super::*;

use crate::tests::read_edits::apply_read_edits;
use crate::tests::read_edits::ReadEdits;
use crate::prelude::CorrectRNAReadsOptions;


// A 200 bp contig with a two-intron gene, 1-based:
//   exon A 1..60, intron 1 61..100, exon B 101..130, intron 2 131..170, exon C 171..200.
//
// Deliberate landmarks: exon A ends in the homopolymer `AAA` at 57..59 (the boundary-contraction
// case); intron 2 and exon C both start with `G` (1 bp of microhomology, the `114=` vs `113=81N`
// overhang case); intron 2's interior differs from exon C everywhere the tests look, so a spliced
// match can never be mistaken for a genomic one.
const EXON_A: &str = "ACGTACGGTTCAGCATTGGAACGTTCGATCCGGATAGCTTAGCGTTACGGATCCATAAAG";
const INTRON_1: &str = "GTGAGTCCTTGCGGACTTAAACCGGTTAACCGGTTAACAG";
const EXON_B: &str = "CATCGCTATCGGATTACCAGTTGACCATGA";
const INTRON_2: &str = "GTAAGGCACGCTTGCAAGGTCTTGCAAGGACTTCCTGCAG";
const EXON_C: &str = "GGCTTACGATTCAGGATCCGTATTGGCTAA";

const CHROMOSOME: u16 = 7;


fn fixture() -> (TempDir, FastaMap) {
    let contig: String = format!("{EXON_A}{INTRON_1}{EXON_B}{INTRON_2}{EXON_C}");
    assert_eq!(contig.len(), 200);
    let directory: TempDir = tempdir().unwrap();
    let fasta_file: String = directory
        .path()
        .join("clip_repair_fixture.fa")
        .to_str()
        .unwrap()
        .to_string();
    std::fs::write(&fasta_file, format!(">chrT\n{contig}\n")).unwrap();
    (directory, FastaMap::new(&fasta_file))
}


fn run(read_start: u32, length: u32, leading: bool, anchor: u32, strand: Strand) -> TerminalSoftClipRun {
    TerminalSoftClipRun {
        read_start,
        length,
        leading,
        chromosome_id: CHROMOSOME,
        anchor,
        strand
    }
}


#[test]
fn test_contig_lengths_come_from_the_fasta_map() {
    let (_directory, fasta_map) = fixture();
    assert_eq!(fasta_map.get_lengths().get("chrT"), Some(&200));
}


/// Strict anchoring is the whole point: a skipped reference base at the boundary is an
/// `ExpectedMissing` that costs — never a free shift. The far end is free, and a minimal path
/// never ends by restoring bases past the read's end.
#[test]
fn test_align_clip_is_anchored_strictly_and_free_only_at_the_far_end() {
    // clip = expected minus its first base: one restore at the boundary, then matches.
    let (cost, ops) = align_clip(b"CGT", b"ACGTAA");
    assert_eq!(cost, 1);
    assert_eq!(
        ops,
        vec![
            ClipAlignmentOp::ExpectedMissing,
            ClipAlignmentOp::Match,
            ClipAlignmentOp::Match,
            ClipAlignmentOp::Match
        ]
    );

    // The continuation being longer than the clip costs nothing.
    let (cost, ops) = align_clip(b"AC", b"ACGTGT");
    assert_eq!(cost, 0);
    assert_eq!(ops, vec![ClipAlignmentOp::Match, ClipAlignmentOp::Match]);
}


/// A clip that exactly matches the genomic continuation needs no repair: the bases are correct
/// and merely unplaceable, so the run resolves to `Trimmable` — kept with trimming off, removed
/// with trimming on — never rewritten.
#[test]
fn test_exact_genomic_continuation_resolves_trimmable() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Aligned across 21..40, clipped on `AGCGTT` = reference 41..46.
    let read: String = format!("{}AGCGTT", &EXON_A[20..40]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 6, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Trimmable);
}


/// A mismatched clip base is substituted with the reference base of the continuation.
#[test]
fn test_genomic_mismatch_is_substituted() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Reference 41..46 is AGCGTT; the clip reads AGTGTT — one mismatch at reference 43.
    let read: String = format!("{}AGTGTT", &EXON_A[20..40]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 6, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 22,
            chromosome_id: CHROMOSOME,
            reference_position: 43,
            quality: 0,
            kind: ReadEditKind::Substitution { reference_base: b'C' }
        }])
    );
}


/// The canonical category-A clip: the read is one base short in the `AAA` run at the boundary,
/// SAM cannot end an alignment with `D`, so the aligner clips the preceding bases — all of them
/// exact reference. The repair restores the missing base at the clip boundary, in read space.
#[test]
fn test_leading_forward_homopolymer_boundary_deletion_is_restored() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Reference 50..60 is G ATCCAT AAA G; the read lost one A and aligns from 58, clipping
    // `GATCCAT` (= reference 50..56, exact). Anchor = alignment start - 1 = 57.
    let read: &str = "GATCCATAAG";
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(0, 7, true, 57, Strand::Forward),
        read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    let expected: Vec<ReadEdit> = vec![ReadEdit {
        read_position: 7,
        chromosome_id: CHROMOSOME,
        reference_position: 57,
        quality: 0,
        kind: ReadEditKind::Deletion { reference_bases: b"A".to_vec() }
    }];
    assert_eq!(resolution, ClipRunResolution::Repair(expected.clone()));

    // Through the emit pass: the restored base completes the homopolymer.
    let alignment: ReadEdits = ReadEdits {
        clip_edits: expected,
        clip_trims: vec![],
        end_trims: vec![],
        ..Default::default()
    };
    let (sequence, quality, num_corrections) = apply_read_edits(
        read.as_bytes(),
        &vec![40u8; read.len()],
        &alignment
    );
    assert_eq!(&sequence, b"GATCCATAAAG");
    assert_eq!(num_corrections, 1);
    assert_eq!(sequence.len(), quality.len());
}


/// The same boundary contraction on the reverse strand: the run sits at the read's trailing end,
/// the walk runs leftward on the reference, and the restored base arrives complemented into read
/// orientation.
#[test]
fn test_trailing_reverse_boundary_deletion_is_restored_complemented() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Same molecule as above, sequenced the other way: read = revcomp(AAG) + revcomp(GATCCAT).
    let read: &str = "CTTATGGATC";
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(3, 7, false, 57, Strand::Reverse),
        read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    let expected: Vec<ReadEdit> = vec![ReadEdit {
        read_position: 3,
        chromosome_id: CHROMOSOME,
        reference_position: 57,
        quality: 0,
        kind: ReadEditKind::Deletion { reference_bases: b"T".to_vec() }
    }];
    assert_eq!(resolution, ClipRunResolution::Repair(expected.clone()));

    let alignment: ReadEdits = ReadEdits {
        clip_edits: expected,
        clip_trims: vec![],
        end_trims: vec![],
        ..Default::default()
    };
    let (sequence, _, _) = apply_read_edits(
        read.as_bytes(),
        &vec![40u8; read.len()],
        &alignment
    );
    // Reverse-complemented back to reference orientation, the homopolymer is whole again.
    assert_eq!(&*reverse_complement(std::str::from_utf8(&sequence).unwrap()), "GATCCATAAAG");
}


/// The category-B clip with 1 bp of microhomology: the alignment overhangs the donor by one base
/// (intron 2 and exon C both start `G`), so the true continuation starts one base into the
/// acceptor exon. A clip that matches it exactly is already correct — `Trimmable`, never
/// "repaired" against the abutting intron.
#[test]
fn test_spliced_overhang_continuation_matches_without_edits() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Leading run on a Reverse record: anchor = alignment end = 131, one base past the donor
    // (intron 2 = 131..170). Expected continuation = reference from 172, complemented:
    // revcomp(GCTTAC) = GTAAGC.
    let read: String = format!("GTAAGC{}", "CCCCC");
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(0, 6, true, 131, Strand::Reverse),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Trimmable);
}


/// The same geometry with one bad base: repaired against the acceptor exon across the intron,
/// never against the intron itself.
#[test]
fn test_spliced_overhang_mismatch_is_repaired_across_the_intron() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // One error in the fragment: read carries GTAATC where the continuation says GTAAGC. In
    // outward order the mismatch sits at reference 173 (exon C), whose read-oriented base is G.
    let read: String = format!("GTAATC{}", "CCCCC");
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(0, 6, true, 131, Strand::Reverse),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 4,
            chromosome_id: CHROMOSOME,
            reference_position: 173,
            quality: 0,
            kind: ReadEditKind::Substitution { reference_base: b'G' }
        }])
    );
}


/// Donor undershoot needs no special case: the walk crosses the remaining exon bases, reaches the
/// donor, and branches across the intron on its own.
#[test]
fn test_spliced_undershoot_walks_to_the_donor_and_repairs_beyond_it() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Alignment ends at 127, three bases short of the donor; the clip spans the exon remainder
    // (TGA = 128..130) plus the start of exon C. Exact clip first:
    let aligned: &str = &EXON_B[0..20]; // read positions 0..19, content irrelevant to the repair
    let read_exact: String = format!("{aligned}TGAGGCTT");
    let exact: ClipRunResolution = resolve_one_clip_run(
        &run(20, 8, false, 127, Strand::Forward),
        &read_exact,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(exact, ClipRunResolution::Trimmable);

    // One error in the acceptor-exon part (reference 173 is C, the read carries A):
    let read_error: String = format!("{aligned}TGAGGATT");
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 8, false, 127, Strand::Forward),
        &read_error,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 25,
            chromosome_id: CHROMOSOME,
            reference_position: 173,
            quality: 0,
            kind: ReadEditKind::Substitution { reference_base: b'C' }
        }])
    );
}


/// A lone terminal artifact base ties drain against substitute; the tie resolves to keeping the
/// base and rewriting it to the continuation — either reading re-aligns clean, and the rule that
/// consumes the most reference is the one the indel-heavy platform prior favors elsewhere.
#[test]
fn test_lone_terminal_artifact_base_is_rewritten_to_the_continuation() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Last aligned base is reference 40 = T; the clip is a second T where the continuation
    // starts A (reference 41).
    let read: &str = "GGGGGGGGGGT";
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(10, 1, false, 40, Strand::Forward),
        read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 10,
            chromosome_id: CHROMOSOME,
            reference_position: 41,
            quality: 0,
            kind: ReadEditKind::Substitution { reference_base: b'A' }
        }])
    );
}


/// A multi-base clip that starts with a duplicated copy of the last aligned base — the measured
/// 5'-terminus artifact — drains the duplicate, because there the drain is strictly cheaper than
/// any substitution reading.
#[test]
fn test_duplicated_base_before_matching_continuation_is_drained() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Last aligned base is reference 40 = T; the clip reads T + AGC, and AGC = reference 41..43.
    let read: &str = "GGGGGGGGGGTAGC";
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(10, 4, false, 40, Strand::Forward),
        read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 10,
            chromosome_id: CHROMOSOME,
            reference_position: 40,
            quality: 0,
            kind: ReadEditKind::Insertion { length: 1 }
        }])
    );
}


/// A boundary contraction across a donor: the clip matches the acceptor exon minus its first `G`
/// (exon C starts `GG`), so one base is restored at the boundary across the intron.
#[test]
fn test_spliced_boundary_deletion_is_restored_across_the_intron() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = format!("{}GCTT", &EXON_B[0..20]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 4, false, 130, Strand::Forward),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 20,
            chromosome_id: CHROMOSOME,
            reference_position: 171,
            quality: 0,
            kind: ReadEditKind::Deletion { reference_bases: b"G".to_vec() }
        }])
    );
}


/// Two continuations at the same cost that would emit different edits are no license to repair.
/// With a second (synthetic) intron sharing the donor, the same clip is one restore against one
/// acceptor and one drain against the other — unrepairable, hence `Trimmable`.
#[test]
fn test_equal_cost_continuations_with_different_edits_resolve_trimmable() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = format!("{}GCTT", &EXON_B[0..20]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 4, false, 130, Strand::Forward),
        &read,
        "chrT",
        &[(131, 170), (131, 160)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Trimmable);
}


/// Low identity is signal, not error, so it is never rewritten — with no called variant behind
/// it, it is `Trimmable`.
#[test]
fn test_low_identity_clip_is_never_rewritten() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = format!("{}GGGGGGCCCCCC", &EXON_A[20..40]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 12, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Trimmable);
    assert_eq!(clip_repair_edit_budget(12), 3);
}


/// A structural-scale run is never repairable. Unlicensed it is `Trimmable`; standing on a
/// called variant it is `Untouched` — the licensing check runs before the length gate, which is
/// what keeps a CALLED breakend's fusion arm safe at any scale.
#[test]
fn test_structural_scale_clip_licensing_precedes_the_length_gate() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = "A".repeat(60);
    let clip_run: TerminalSoftClipRun = run(5, 51, false, 40, Strand::Forward);

    let unlicensed: ClipRunResolution = resolve_one_clip_run(
        &clip_run,
        &read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(unlicensed, ClipRunResolution::Trimmable);

    let protected: HashSet<(u16, u32)> = HashSet::from([(CHROMOSOME, 40)]);
    let licensed: ClipRunResolution = resolve_one_clip_run(
        &clip_run,
        &read,
        "chrT",
        &[],
        Some(&protected),
        &fasta_map,
        &no_trim
    );
    assert_eq!(licensed, ClipRunResolution::Untouched);
}


/// An indel in a repeat has no canonical placement: the caller may left-align the deletion in
/// the `AAA` run (protecting 56..58) while the boundary-anchored repair spells the same event at
/// 59. Exact-position matching would wave that repair through and erase the called variant's
/// evidence; the proximity window leaves the run untouched.
#[test]
fn test_equivalently_placed_indel_near_a_called_variant_is_untouched() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Read carries the deletion the caller called: two A's aligned at 57..58, the trailing G
    // clipped. Unprotected, the repair restores the A — spelled at 59, the run's right end.
    let read: &str = "GGGGGGGGGGG";
    let unprotected: ClipRunResolution = resolve_one_clip_run(
        &run(10, 1, false, 58, Strand::Forward),
        read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        unprotected,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 10,
            chromosome_id: CHROMOSOME,
            reference_position: 59,
            quality: 0,
            kind: ReadEditKind::Deletion { reference_bases: b"A".to_vec() }
        }])
    );

    // The caller's left-aligned spelling protects 56..58 — none of which the repair touches
    // exactly. The window catches it anyway, and the run is neither repaired nor trimmed.
    let protected: HashSet<(u16, u32)> =
        HashSet::from([(CHROMOSOME, 56), (CHROMOSOME, 57), (CHROMOSOME, 58)]);
    let vetoed: ClipRunResolution = resolve_one_clip_run(
        &run(10, 1, false, 58, Strand::Forward),
        read,
        "chrT",
        &[],
        Some(&protected),
        &fasta_map,
        &no_trim
    );
    assert_eq!(vetoed, ClipRunResolution::Untouched);
}


/// A called breakend protects only its two anchors, and a fusion-arm clip stands exactly on one
/// of them. The anchor is licensed up front, so the run is `Untouched` — not repaired, and with
/// trimming enabled, not trimmed either.
#[test]
fn test_clip_anchored_on_a_called_breakend_is_untouched() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = format!("{}AGTGTT", &EXON_A[20..40]);
    let protected: HashSet<(u16, u32)> = HashSet::from([(CHROMOSOME, 40)]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 6, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        Some(&protected),
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Untouched);
}


/// A repair resting exactly on a called variant's position is likewise left untouched.
#[test]
fn test_repair_touching_a_protected_position_is_untouched() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    let read: String = format!("{}AGTGTT", &EXON_A[20..40]);
    let protected: HashSet<(u16, u32)> = HashSet::from([(CHROMOSOME, 43)]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 6, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        Some(&protected),
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Untouched);
}


/// A maximal in-budget boundary deletion must be expressible: the continuation padding exceeds
/// the largest edit budget, so a 9 bp boundary deletion on a 40 bp clip is restored whole rather
/// than being mis-scored against a truncated continuation (which would drain a real base).
#[test]
fn test_large_boundary_deletion_is_restored_whole() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Clip = reference 50..89 exactly; reference 41..49 (AGCGTTACG) is missing from the read.
    let contig: String = format!("{EXON_A}{INTRON_1}{EXON_B}{INTRON_2}{EXON_C}");
    let clip: &str = &contig[49..89];
    let read: String = format!("{}{}", &EXON_A[20..40], clip);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 40, false, 40, Strand::Forward),
        &read,
        "chrT",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 20,
            chromosome_id: CHROMOSOME,
            reference_position: 41,
            quality: 0,
            kind: ReadEditKind::Deletion { reference_bases: b"AGCGTTACG".to_vec() }
        }])
    );
}


/// A continuation running into an assembly gap must not write `N` into the read: an in-budget
/// alignment whose repairs would emit non-ACGT bases resolves `Trimmable`, never `Repair`.
#[test]
fn test_repair_never_writes_non_acgt_bases() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let directory: TempDir = tempdir().unwrap();
    let fasta_file: String = directory
        .path()
        .join("gap_fixture.fa")
        .to_str()
        .unwrap()
        .to_string();
    // 20 real bases, then a 20 bp assembly gap.
    std::fs::write(&fasta_file, format!(">chrU\n{}{}\n", &EXON_A[0..20], "N".repeat(20))).unwrap();
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);

    // Clip = reference 9..20 exactly, then four read bases standing over the gap: 12 matches
    // plus 4 N-mismatches is cost 4 = the budget for a 16 bp clip, so only the alphabet gate
    // stands between those bases and literal Ns at Q30.
    let read: String = format!("GGGGGGGG{}TTTT", &EXON_A[8..20]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(8, 16, false, 8, Strand::Forward),
        &read,
        "chrU",
        &[],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(resolution, ClipRunResolution::Trimmable);
}


/// Clip edits go through the emit pass with the fixed repair quality planning stamped on them;
/// no site evidence is looked up for a clip anchor.
#[test]
fn test_correct_one_read_applies_clip_edits_with_fixed_quality() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let clip_edits: Vec<ReadEdit> = vec![
        ReadEdit {
            read_position: 4,
            chromosome_id: 0,
            reference_position: 105,
            quality: options.soft_clip_repair_base_quality,
            kind: ReadEditKind::Substitution { reference_base: b'A' }
        },
        ReadEdit {
            read_position: 8,
            chromosome_id: 0,
            reference_position: 200,
            quality: options.soft_clip_repair_base_quality,
            kind: ReadEditKind::Deletion { reference_bases: b"T".to_vec() }
        },
    ];
    let alignment: ReadEdits = ReadEdits {
        clip_edits: clip_edits.clone(),
        ..Default::default()
    };

    let (sequence, quality, num_corrections) = apply_read_edits(
        b"ACGTCCGT",
        &[40u8; 8],
        &alignment
    );
    assert_eq!(&sequence, b"ACGTACGTT");
    assert_eq!(num_corrections, 2);
    assert_eq!(sequence.len(), quality.len());
    assert_eq!(quality[4], b'!' + options.soft_clip_repair_base_quality);
    assert_eq!(quality[8], b'!' + options.soft_clip_repair_base_quality);
}


/// A repairable clip standing entirely beyond a junction is a next-exon fragment: with trimming
/// off it is repaired whole (shipped behavior), but with trimming on the whole run is trimmed —
/// the repaired fragment would sit below the aligner's placeable size and come straight back as a
/// clip. This is the `chunk_0000/918`-family case that motivated the split.
#[test]
fn test_repairable_post_junction_fragment_is_trimmed_whole_when_enabled() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Same geometry as the spliced-overhang repair test: one error in a 6 bp fragment whose only
    // in-budget continuation starts across intron 2 (overhang start, a jump before base one).
    let read: String = format!("GTAATC{}", "CCCCC");
    let repaired: ClipRunResolution = resolve_one_clip_run(
        &run(0, 6, true, 131, Strand::Reverse),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert!(matches!(repaired, ClipRunResolution::Repair(_)));

    let trimmed: ClipRunResolution = resolve_one_clip_run(
        &run(0, 6, true, 131, Strand::Reverse),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &options
    );
    assert_eq!(
        trimmed,
        ClipRunResolution::RepairAndTrim {
            edits: vec![],
            trim: ReadEdit {
                read_position: 0,
                chromosome_id: CHROMOSOME,
                reference_position: 131,
                quality: 0,
                kind: ReadEditKind::Insertion { length: 6 }
            }
        }
    );
}


/// A clip that continues the current exon and then crosses the junction is split at the junction:
/// the genomic prefix is repaired and kept (it extends the alignment flush on re-alignment), the
/// short post-junction remainder is trimmed. Trimming off keeps the whole repair, as shipped.
#[test]
fn test_mid_clip_junction_split_keeps_the_prefix_and_trims_the_remainder() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let no_trim: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_unsupported_soft_clips: false, ..options.clone() };
    let (_directory, fasta_map) = fixture();
    // Alignment ends at 127; the clip spans the exon remainder TGA (128..130, read carries TCA —
    // one error) and then 5 bases of exon C across intron 2.
    let read: String = format!("{}TCAGGCTT", &EXON_B[0..20]);
    let prefix_repair: ReadEdit = ReadEdit {
        read_position: 21,
        chromosome_id: CHROMOSOME,
        reference_position: 129,
        quality: 0,
        kind: ReadEditKind::Substitution { reference_base: b'G' }
    };

    let whole: ClipRunResolution = resolve_one_clip_run(
        &run(20, 8, false, 127, Strand::Forward),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &no_trim
    );
    assert_eq!(whole, ClipRunResolution::Repair(vec![prefix_repair.clone()]));

    let split: ClipRunResolution = resolve_one_clip_run(
        &run(20, 8, false, 127, Strand::Forward),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &options
    );
    assert_eq!(
        split,
        ClipRunResolution::RepairAndTrim {
            edits: vec![prefix_repair],
            trim: ReadEdit {
                read_position: 23,
                chromosome_id: CHROMOSOME,
                reference_position: 127,
                quality: 0,
                kind: ReadEditKind::Insertion { length: 5 }
            }
        }
    );
}


/// A post-junction remainder at or above the placeable size is not split: the aligner will pay
/// the intron for a terminal block that long, so the whole repair survives re-alignment.
#[test]
fn test_placeable_post_junction_remainder_is_kept_whole() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let (_directory, fasta_map) = fixture();
    // Same prefix error as above, but the clip carries 16 bases of exon C — placeable.
    let read: String = format!("{}TCA{}", &EXON_B[0..20], &EXON_C[0..16]);
    let resolution: ClipRunResolution = resolve_one_clip_run(
        &run(20, 19, false, 127, Strand::Forward),
        &read,
        "chrT",
        &[(131, 170)],
        None,
        &fasta_map,
        &options
    );
    assert_eq!(
        resolution,
        ClipRunResolution::Repair(vec![ReadEdit {
            read_position: 21,
            chromosome_id: CHROMOSOME,
            reference_position: 129,
            quality: 0,
            kind: ReadEditKind::Substitution { reference_base: b'G' }
        }])
    );
}


/// A clip trim drains its whole run through the emit pass — deliberately without the
/// structural-scale cap, since trimming an uncalled breakend's arm is the point — and the
/// quality string stays in step.
#[test]
fn test_correct_one_read_applies_clip_trims_beyond_the_scale_cap() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let trim_length: usize = options.max_correctable_event_len + 1;
    let mut sequence: Vec<u8> = b"ACGTACGT".to_vec();
    sequence.extend(std::iter::repeat(b'A').take(trim_length));
    let quality: Vec<u8> = vec![40u8; sequence.len()];

    let alignment: ReadEdits = ReadEdits {
        clip_edits: vec![],
        clip_trims: vec![ReadEdit {
            read_position: 8,
            chromosome_id: 0,
            reference_position: 300,
            quality: 0,
            kind: ReadEditKind::Insertion { length: trim_length }
        }],
        end_trims: vec![],
        ..Default::default()
    };

    let (trimmed, trimmed_quality, num_corrections) = apply_read_edits(
        &sequence,
        &quality,
        &alignment
    );
    assert_eq!(&trimmed, b"ACGTACGT", "the trim did not remove the whole clip run");
    assert_eq!(trimmed.len(), trimmed_quality.len());
    assert_eq!(num_corrections, 1);
}
