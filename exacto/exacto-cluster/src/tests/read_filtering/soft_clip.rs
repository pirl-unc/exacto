use std::fs;
use std::path::Path;

use super::*;


/// A terminal soft clip is the end of a degraded read when it spells the reference across an
/// intron next to it.
///
/// The reads sit on TP53 (chr17, reverse strand) around the intron 7,674,972-7,675,052. The
/// exon above it starts with CATCGCTATC at 7,675,053 and the exon below it ends with
/// GGGCCAGAC at 7,674,971. The intron starts with C as the exon above does, so the aligner
/// runs a read that stops inside the exon above one base into the intron, to 7,674,972, and
/// clips the rest.
///
///   Clip along the reference   Hangs off    Read end   Identified
///   ATCGCTATC                  7,674,972    5'         yes
///   ATCGCTA                    7,674,972    5'         yes
///   ATCGaTATC                  7,674,972    5'         yes, 1 edit in 9 bases
///   ATCGaTA                    7,674,972    5'         no, no edit in 7 bases
///   GGGTTTAAA                  7,674,972    5'         no
///   ATCGCTATC                  7,674,972    inside     no, an insertion within the read
///   GGGCCAGAC                  7,675,053    3'         yes
///   ATCGCTATC                  7,674,872    5'         no, no intron within 3 bases
///
/// That is with a boundary within 3 bases of the end of the alignment and 1 edit for every 8
/// bases of the clip. The clip at 7,674,972 hangs off one base from the boundary, 7,674,971:
///
///   Clip along the reference   Boundary within   Bases per edit   Identified
///   ATCGCTATC                  1                 8                yes
///   ATCGCTATC                  0                 8                no, the boundary is 1 base away
///   ATCGaTATC                  3                 0                no, no edit is allowed
///   ATCGCTATC                  3                 0                yes, it has no edit
///   ATCGaTA                    3                 4                yes, 1 edit in 7 bases
#[test]
fn is_soft_clip_across_intron_returns_matches() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([(Box::from("chr17"), 0u16)]);

    let summaries: Vec<RNAReadCharacterizationSummary> = vec![
        RNAReadCharacterizationSummary {
            read_id: 0,
            chromosome: Some(0),
            reference_start: 7_675_236,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_053, 7_675_236), (0, 7_674_859, 7_674_971)],
            offset: 0,
            length: 0
        }
    ];
    let intron_boundary_index: HashMap<(u16, u32), Vec<(u32, u32)>> = build_intron_boundary_index(&summaries);
    assert_eq!(
        intron_boundary_index,
        HashMap::from([
            ((0, 7_674_972), vec![(7_674_972, 7_675_052)]),
            ((0, 7_675_052), vec![(7_674_972, 7_675_052)])
        ])
    );

    // (read position 1, read position 2, position 1, position 2, sequence on the read, identified)
    let read_length: u32 = 1_000;
    let cases: Vec<(u32, u32, u32, u32, &str, bool)> = vec![
        (0, 8, 7_674_972, 7_674_973, "GATAGCGAT", true),
        (0, 6, 7_674_972, 7_674_973, "TAGCGAT", true),
        (0, 8, 7_674_972, 7_674_973, "GATATCGAT", true),
        (0, 6, 7_674_972, 7_674_973, "TATCGAT", false),
        (0, 8, 7_674_972, 7_674_973, "TTTAAACCC", false),
        (400, 408, 7_674_972, 7_674_973, "GATAGCGAT", false),
        (991, 999, 7_675_052, 7_675_053, "GTCTGGCCC", true),
        (0, 8, 7_674_872, 7_674_873, "GATAGCGAT", false)
    ];
    for (read_position_1, read_position_2, position_1, position_2, sequence, expected) in cases {
        let variant_record: VariantRecord = VariantRecord::new(
            0,
            read_position_1,
            read_position_2,
            GraphOperation::new(
                0,
                position_1,
                Strand::Reverse,
                GraphOperationType::Downstream,
                0,
                position_2,
                Strand::Reverse,
                GraphOperationType::Upstream,
                sequence.into(),
                VariantType::Insertion
            )
        );
        assert_eq!(
            is_soft_clip_across_intron(
                &variant_record,
                read_length,
                &[],
                &intron_boundary_index,
                &chromosome_names_map,
                &fasta_map,
                3,
                8
            ),
            expected,
            "{} at {}-{}, read positions {}-{}",
            sequence,
            position_1,
            position_2,
            read_position_1,
            read_position_2
        );
    }

    // (sequence on the read, max boundary distance, number of bases per edit, identified)
    let cases: Vec<(&str, u32, u32, bool)> = vec![
        ("GATAGCGAT", 1, 8, true),
        ("GATAGCGAT", 0, 8, false),
        ("GATATCGAT", 3, 0, false),
        ("GATAGCGAT", 3, 0, true),
        ("TATCGAT", 3, 4, true)
    ];
    for (sequence, max_boundary_distance, num_bases_per_edit, expected) in cases {
        let variant_record: VariantRecord = VariantRecord::new(
            0,
            0,
            sequence.len() as u32 - 1,
            GraphOperation::new(
                0,
                7_674_972,
                Strand::Reverse,
                GraphOperationType::Downstream,
                0,
                7_674_973,
                Strand::Reverse,
                GraphOperationType::Upstream,
                sequence.into(),
                VariantType::Insertion
            )
        );
        assert_eq!(
            is_soft_clip_across_intron(
                &variant_record,
                read_length,
                &[],
                &intron_boundary_index,
                &chromosome_names_map,
                &fasta_map,
                max_boundary_distance,
                num_bases_per_edit
            ),
            expected,
            "{} with a boundary within {} and {} bases per edit",
            sequence,
            max_boundary_distance,
            num_bases_per_edit
        );
    }
}


/// A terminal soft clip is the end of a degraded read when it spells the start of an exon the
/// read itself holds.
///
/// The read is of a TP53 transcript (chr17, reverse strand) with a tandem duplication of four
/// exons, from 7,673,701-7,673,837 to 7,675,053-7,675,236. After the last of them the read
/// goes round to the first again. It starts 68 bases into that copy, too few for the aligner to
/// place, so the alignment ends at the intron 7,675,237-7,675,993 and the rest is clipped. The
/// intron starts with TC as the exon at 7,673,701 does, so the aligner runs the read two bases
/// into it, to 7,675,238, and the clip starts at 7,673,703.
///
///   Clip along the reference   Introns the read splices   Identified
///   7,673,703-7,673,768        7,673,609-7,673,700        yes
///   7,673,703-7,673,722        7,673,609-7,673,700        yes
///   7,673,703-7,673,768        7,674,291-7,674,858        no, the exon above it starts otherwise
///   7,673,703-7,673,768        none                       no
#[test]
fn is_soft_clip_across_intron_returns_matches_for_exons_of_read() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([(Box::from("chr17"), 0u16)]);

    let summaries: Vec<RNAReadCharacterizationSummary> = vec![
        RNAReadCharacterizationSummary {
            read_id: 0,
            chromosome: Some(0),
            reference_start: 7_676_128,
            reference_end: 7_675_053,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_993, 7_675_237, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_994, 7_676_128), (0, 7_675_053, 7_675_236)],
            offset: 0,
            length: 0
        }
    ];
    let intron_boundary_index: HashMap<(u16, u32), Vec<(u32, u32)>> = build_intron_boundary_index(&summaries);

    // (sequence on the read, introns the read splices, identified)
    let cases: Vec<(&str, Vec<(u16, u32, u32)>, bool)> = vec![
        (
            "AGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCG",
            vec![(0, 7_673_609, 7_673_700)],
            true
        ),
        ("CCCCCAGGGAGCACTAAGCG", vec![(0, 7_673_609, 7_673_700)], true),
        (
            "AGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCG",
            vec![(0, 7_674_291, 7_674_858)],
            false
        ),
        (
            "AGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCG",
            Vec::new(),
            false
        )
    ];
    for (sequence, read_introns, expected) in cases {
        let variant_record: VariantRecord = VariantRecord::new(
            0,
            0,
            sequence.len() as u32 - 1,
            GraphOperation::new(
                0,
                7_675_238,
                Strand::Reverse,
                GraphOperationType::Downstream,
                0,
                7_675_239,
                Strand::Reverse,
                GraphOperationType::Upstream,
                sequence.into(),
                VariantType::Insertion
            )
        );
        assert_eq!(
            is_soft_clip_across_intron(
                &variant_record,
                2_063,
                &read_introns,
                &intron_boundary_index,
                &chromosome_names_map,
                &fasta_map,
                3,
                8
            ),
            expected,
            "{} with read introns {:?}",
            sequence,
            read_introns
        );
    }
}


#[test]
fn is_within_edits_of_a_prefix_returns_matches_of_edit_distance_per_prefix() {
    let mut state: u64 = 42;
    let mut next = || -> usize {
        state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) as usize
    };
    let mut num_within: usize = 0;
    for trial in 0..3_000usize {
        let clip: Vec<u8> = (0..1 + next() % 40).map(|_| b"ACGT"[next() % 4]).collect();
        let mut target: Vec<u8> = if trial % 3 == 0 {
            (0..next() % 50).map(|_| b"ACGT"[next() % 4]).collect()
        } else {
            clip.clone()
        };
        if trial % 3 != 0 {
            for _ in 0..next() % 7 {
                let position: usize = next() % (target.len() + 1);
                match next() % 3 {
                    0 if position < target.len() => target[position] = b"ACGT"[next() % 4],
                    1 if position < target.len() => { target.remove(position); },
                    _ => target.insert(position, b"ACGT"[next() % 4])
                }
            }
            target.extend((0..next() % 10).map(|_| b"ACGT"[next() % 4]));
        }
        let max_edits: usize = (next() % 9).min(clip.len());
        let min_length: usize = clip.len() - max_edits;
        let max_length: usize = (clip.len() + max_edits).min(target.len());
        let clip_text: &str = std::str::from_utf8(&clip).unwrap();
        let target_text: &str = std::str::from_utf8(&target).unwrap();
        let expected: bool = (min_length..=max_length)
            .any(|length| edit_distance::edit_distance(clip_text, &target_text[..length]) <= max_edits);
        num_within += usize::from(expected);
        assert_eq!(
            is_within_edits_of_a_prefix(&clip, &target, min_length, max_length, max_edits),
            expected,
            "clip {clip_text}, target {target_text}, {max_edits} edits"
        );
    }
    assert!(num_within > 500 && num_within < 2_500, "{num_within} of 3000 within");
}
