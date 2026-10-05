use std::collections::{HashMap, HashSet};
use std::path::Path;

use exacto_cluster::prelude::{RNAReadCluster, RNAReadClusterSet};
use tempfile::tempdir;

use super::*;

use crate::prelude::{correct_rna_reads, CorrectRNAReadsOptions};


// The clip-repair landmarks, reused: exon A 1..60, exon B 101..130, exon C 171..200 on a
// 200 bp contig, so intron 1 is 61..100 and intron 2 is 131..170.


/// `exons` in transcript order, exon 1 first, so a minus-strand transcript lists them high
/// to low.
fn transcript(strand: Strand, exons: &[(u32, u32)]) -> Transcript {
    let start: u32 = exons.iter().map(|exon| exon.0).min().unwrap();
    let end: u32 = exons.iter().map(|exon| exon.1).max().unwrap();
    let mut transcript: Transcript = Transcript::new(
        "GENE", "TX", "test", "chrT", start, end, strand.clone(), 1, "TX", "protein_coding", "1",
        HashSet::new()
    );
    for (index, (exon_start, exon_end)) in exons.iter().enumerate() {
        transcript.add_exon(Exon::new(
            "GENE", "TX", &format!("EX{}", index + 1), "test", "chrT",
            *exon_start, *exon_end, strand.clone(), 1, (index + 1) as u16
        ));
    }
    transcript
}


fn block(start: u32, end: u32, strand: Strand, read_start: u32, read_end: u32) -> TranscriptModelExon {
    TranscriptModelExon::new(0, start, end, strand, 0, read_start, read_end)
}


#[test]
fn test_tail_past_an_internal_exon_is_cut_at_the_exon_end() {
    let blocks: Vec<TranscriptModelExon> = vec![
        block(21, 60, Strand::Forward, 0, 39),
        block(101, 145, Strand::Forward, 40, 84)
    ];
    let terminals = terminal_blocks(&blocks);
    assert_eq!(terminals.len(), 2);
    let tail = terminals.iter().find(|terminal| terminal.end == ReadEnd::Tail).unwrap();
    assert_eq!(tail.intron, (61, 100));
    assert!(tail.rightward);

    let three_exon: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 130), (171, 200)]);
    assert_eq!(judge_end(tail, &three_exon), Some(EndVerdict::Boundary(130)));
    assert_eq!(resolve_cut_boundary(tail, &[&three_exon]), Some(130));

    // Read 40 is placed at 101, so 131 is read 70: fifteen bases go.
    let placements: Vec<(u32, u32)> = (40..=84).map(|position| (position, position + 61)).collect();
    assert_eq!(cut_from_placements(&placements, 130, true, ReadEnd::Tail, 85), Some((70, 15)));

    // The head sits in exon 1, a transcript end, and is never cut.
    let head = terminals.iter().find(|terminal| terminal.end == ReadEnd::Head).unwrap();
    assert_eq!(judge_end(head, &three_exon), Some(EndVerdict::TranscriptEnd));
    assert_eq!(resolve_cut_boundary(head, &[&three_exon]), None);
}


#[test]
fn test_reverse_strand_tail_is_cut_at_the_exon_start() {
    // A sense read on a minus-strand gene: the head sits at the high coordinates.
    let blocks: Vec<TranscriptModelExon> = vec![
        block(171, 200, Strand::Reverse, 0, 29),
        block(86, 130, Strand::Reverse, 30, 74)
    ];
    let tail = terminal_blocks(&blocks)
        .into_iter()
        .find(|terminal| terminal.end == ReadEnd::Tail)
        .unwrap();
    assert_eq!(tail.intron, (131, 170));
    assert!(!tail.rightward);

    let minus: Transcript = transcript(Strand::Reverse, &[(171, 200), (101, 130), (1, 60)]);
    assert_eq!(resolve_cut_boundary(&tail, &[&minus]), Some(101));
    // A plus-strand transcript with the same introns says nothing about a reverse-strand read.
    let plus: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 130), (171, 200)]);
    assert_eq!(judge_end(&tail, &plus), None);

    // Read 30 is placed at 130 and read 59 at 101, so read 60 is the first base below it.
    let placements: Vec<(u32, u32)> = (30..=74).map(|position| (position, 160 - position)).collect();
    assert_eq!(cut_from_placements(&placements, 101, false, ReadEnd::Tail, 75), Some((60, 15)));
}


/// The mirror of the forward tail: a reverse-strand read's head sits at the high coordinates
/// and leaves its block rightward, so it is cut back to the exon's end, from the read's start.
#[test]
fn test_reverse_strand_head_is_cut_at_the_exon_end() {
    let blocks: Vec<TranscriptModelExon> = vec![
        block(101, 145, Strand::Reverse, 0, 44),
        block(21, 60, Strand::Reverse, 45, 84)
    ];
    let head = terminal_blocks(&blocks)
        .into_iter()
        .find(|terminal| terminal.end == ReadEnd::Head)
        .unwrap();
    assert_eq!(head.intron, (61, 100));
    assert!(head.rightward);

    let minus: Transcript = transcript(Strand::Reverse, &[(171, 200), (101, 130), (1, 60)]);
    assert_eq!(judge_end(&head, &minus), Some(EndVerdict::Boundary(130)));
    assert_eq!(resolve_cut_boundary(&head, &[&minus]), Some(130));

    // Read 0 is placed at 145 and read 14 at 131: the first fifteen bases go.
    let placements: Vec<(u32, u32)> = (0..=44).map(|position| (position, 145 - position)).collect();
    assert_eq!(cut_from_placements(&placements, 130, true, ReadEnd::Head, 85), Some((0, 15)));

    // A block flush with the boundary has nothing beyond it.
    let flush: Vec<(u32, u32)> = (0..=29).map(|position| (position, 130 - position)).collect();
    assert_eq!(cut_from_placements(&flush, 130, true, ReadEnd::Head, 70), None);
}


#[test]
fn test_head_past_an_internal_exon_keeps_an_insertion_anchored_inside_it() {
    let blocks: Vec<TranscriptModelExon> = vec![
        block(90, 130, Strand::Forward, 0, 40),
        block(171, 200, Strand::Forward, 41, 70)
    ];
    let head = terminal_blocks(&blocks)
        .into_iter()
        .find(|terminal| terminal.end == ReadEnd::Head)
        .unwrap();
    assert_eq!(head.intron, (131, 170));
    assert!(!head.rightward);

    let three_exon: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 130), (171, 200)]);
    assert_eq!(resolve_cut_boundary(&head, &[&three_exon]), Some(101));

    // Bases 0..=4 sit below 101; the two bases inserted at anchor 101 stay.
    let placements: Vec<(u32, u32)> = vec![
        (0, 96), (1, 97), (2, 98), (3, 99), (4, 100), (5, 101), (6, 101), (7, 101), (8, 102)
    ];
    assert_eq!(cut_from_placements(&placements, 101, false, ReadEnd::Head, 9), Some((0, 5)));
}


#[test]
fn test_an_alternative_donor_and_a_transcript_end_both_keep_the_read() {
    let blocks: Vec<TranscriptModelExon> = vec![
        block(21, 60, Strand::Forward, 0, 39),
        block(101, 135, Strand::Forward, 40, 74)
    ];
    let tail = terminal_blocks(&blocks)
        .into_iter()
        .find(|terminal| terminal.end == ReadEnd::Tail)
        .unwrap();
    let short_exon: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 130), (171, 200)]);
    let long_exon: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 140), (171, 200)]);
    let two_exon: Transcript = transcript(Strand::Forward, &[(1, 60), (101, 130)]);

    assert_eq!(resolve_cut_boundary(&tail, &[&short_exon]), Some(130));
    // The block stops inside the longer annotated exon.
    assert_eq!(resolve_cut_boundary(&tail, &[&short_exon, &long_exon]), None);
    // One transcript ending at the exon outranks every other verdict.
    assert_eq!(resolve_cut_boundary(&tail, &[&short_exon, &two_exon]), None);
    // No transcript carries the intron.
    assert_eq!(resolve_cut_boundary(&tail, &[]), None);
}


#[test]
fn test_unanchored_ends_are_never_judged() {
    // One block: nothing anchors either end.
    assert!(terminal_blocks(&[block(21, 60, Strand::Forward, 0, 39)]).is_empty());
    // A fusion: the tail's neighbour is on another chromosome.
    let mut other: TranscriptModelExon = block(101, 145, Strand::Forward, 40, 84);
    other.reference_chromosome_id = 1;
    assert!(terminal_blocks(&[block(21, 60, Strand::Forward, 0, 39), other]).is_empty());
    // A back-splice: the tail lies below its neighbour on the forward strand.
    let back_spliced: Vec<TranscriptModelExon> = vec![
        block(101, 130, Strand::Forward, 0, 29),
        block(21, 60, Strand::Forward, 30, 69)
    ];
    assert!(terminal_blocks(&back_spliced).is_empty());
}


/// A trusted operation vetoes a cut when it names a position inside the removed span, or when
/// it is a small variant whose reference span covers it. A structural operation counts only by
/// its anchors: its span may be a whole locus, and filling that in would disable trimming
/// across it.
#[test]
fn test_protection_counts_small_variant_spans_but_only_structural_anchors() {
    let removed: RangeInclusive<u32> = 121..=123;
    let small = |left: u32, right: u32, sequence: &str, kind: VariantType| -> GraphOperation {
        GraphOperation::new(
            0, left, Strand::Forward, GraphOperationType::Downstream,
            0, right, Strand::Forward, GraphOperationType::Upstream, sequence.into(), kind
        )
    };

    // An anchor inside the span.
    assert!(overlaps_protected_operation(&small(121, 123, "C", VariantType::SingleNucleotideVariant), 0, &removed));
    // Anchors outside, span covering: counts for a small variant ...
    assert!(overlaps_protected_operation(&small(120, 124, "CCC", VariantType::MultiNucleotideVariant), 0, &removed));
    assert!(overlaps_protected_operation(&small(110, 130, "", VariantType::Deletion), 0, &removed));
    // ... but not for a structural operation.
    assert!(!overlaps_protected_operation(&GraphOperation::new(
        0, 110, Strand::Forward, GraphOperationType::Downstream,
        0, 130, Strand::Forward, GraphOperationType::Downstream, "".into(), VariantType::Breakpoint
    ), 0, &removed));
    // A structural anchor inside the span counts, whatever the other end is.
    assert!(overlaps_protected_operation(&GraphOperation::new(
        0, 122, Strand::Forward, GraphOperationType::Downstream,
        1, 200, Strand::Forward, GraphOperationType::Upstream, "".into(), VariantType::FusionGene
    ), 0, &removed));
    // Another chromosome, or a small variant entirely outside the span, does not.
    assert!(!overlaps_protected_operation(&GraphOperation::new(
        1, 121, Strand::Forward, GraphOperationType::Downstream,
        1, 123, Strand::Forward, GraphOperationType::Upstream, "C".into(), VariantType::SingleNucleotideVariant
    ), 0, &removed));
    assert!(!overlaps_protected_operation(&small(99, 101, "C", VariantType::SingleNucleotideVariant), 0, &removed));
}


// The fixture in `$EXACTO_TEST_DATA/exacto/exacto-correction/`: the same
// landmarks as above on a 240 bp contig, one three-exon transcript over exons A, B and C, and
// three error-free reads whose ends are the cases the rule must separate.
const EXON_A: &str = "ACGTACGGTTCAGCATTGGAACGTTCGATCCGGATAGCTTAGCGTTACGGATCCATAAAG";
const INTRON_1: &str = "GTGAGTCCTTGCGGACTTAAACCGGTTAACCGGTTAACAG";
const EXON_B: &str = "CATCGCTATCGGATTACCAGTTGACCATGA";
const INTRON_2: &str = "GTAAGGCACGCTTGCAAGGTCTTGCAAGGACTTCCTGCAG";
const EXON_C: &str = "GGCTTACGATTCAGGATCCGTATTGGCTAA";
const DOWNSTREAM: &str = "CGCGTTGCCGGTTCGCCGGATCGCGTTGGCCGATCGGCCG";


/// A gap between two aligned blocks licenses a cut only when a splice event spans it. The
/// tail-overhang read aligned as two records with no `~` between them, a same-chromosome split
/// alignment, yields the same two blocks and the same overhang, and no cut.
#[test]
fn test_same_chromosome_split_alignments_do_not_license_trimming() {
    use noodles_bam as bam;
    use noodles_sam as sam;
    use sam::alignment::io::Write;
    use std::sync::Arc;

    let sequence: String = format!("{}{}{}", &EXON_A[20..60], EXON_B, &INTRON_2[..15]);
    let quality: String = "I".repeat(sequence.len());
    let header: &str = "@HD\tVN:1.6\n@SQ\tSN:chrT\tLN:240\n";
    let spliced: String = format!(
        "{header}tail_overhang\t0\tchrT\t21\t60\t40M40N45M\t*\t0\t0\t{sequence}\t{quality}\tcs:Z::40~gt40ag:45\n"
    );
    let split: String = format!(
        "{header}tail_overhang\t0\tchrT\t21\t60\t40M45S\t*\t0\t0\t{sequence}\t{quality}\tcs:Z::40\n\
         tail_overhang\t2048\tchrT\t101\t60\t40S45M\t*\t0\t0\t{sequence}\t{quality}\tcs:Z::45\n"
    );
    let annotation_file = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-correction/gene_annotations.tsv");
    let annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "test", "v1");
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([("chrT".into(), 0u16)]);

    for (sam_text, expected) in [(spliced, vec![70..85]), (split, Vec::new())] {
        let mut reader = sam::io::Reader::new(sam_text.as_bytes());
        let header = reader.read_header().unwrap();
        let mut writer = bam::io::Writer::new(Vec::new());
        for record in reader.record_bufs(&header) {
            writer.write_alignment_record(&header, &record.unwrap()).unwrap();
        }
        writer.try_finish().unwrap();
        let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
        let mut records: Vec<bam::Record> = Vec::new();
        loop {
            let mut record = bam::Record::default();
            if reader.read_record(&mut record).unwrap() == 0 {
                break;
            }
            records.push(record);
        }
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let model: AlignmentModel = AlignmentModel::new(7, &read_sequence, &base_quality_scores, &records);
        let blocks: Vec<TranscriptModelExon> = identify_transcript_model_exons(&model);
        assert_eq!(blocks.len(), 2, "both alignments must yield the same two blocks");

        let trims: Vec<Range<u32>> = identify_transcript_end_trims(
            &model, &blocks, &chromosome_names_map, &annotator, &HashSet::new()
        );
        assert_eq!(trims, expected);
    }
}


/// End to end over the synthetic BAM: a tail that runs 15 bases past exon B is cut back to it,
/// a head that starts 11 bases before exon B is cut forward to it, and a tail that runs past
/// exon C, the transcript's last exon, is left alone. A called variant inside a removed span
/// keeps that end, and is then applied.
#[test]
fn test_transcript_end_trimming_end_to_end() {
    let directory = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-correction");
    let canonical = |name: &str| -> String {
        std::fs::canonicalize(directory.join(name)).unwrap().to_str().unwrap().to_string()
    };
    let (bam_file, annotation_file) = (canonical("reads.bam"), canonical("gene_annotations.tsv"));
    let annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(&annotation_file, "test", "v1");

    let read_names: HashSet<Box<str>> = ["tail_overhang", "head_overhang", "last_exon_overhang"]
        .into_iter()
        .map(Into::into)
        .collect();
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrT".into(), 0u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.iter().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }
    // Each variant as a one-record call: with a single record the consensus is that record's
    // own operation, so the alignment scores never come into play.
    let build = |variants: Vec<GraphOperation>| -> RNAReadClusterSet {
        let variant_calls: Vec<VariantCall> = variants
            .into_iter()
            .enumerate()
            .map(|(id, variant)| VariantCall::from_variant_records(
                id, HashSet::from([VariantRecord::new(0, 0, 0, variant)]), 2, -4, 4, 2
            ))
            .collect();
        let mut cluster_set = RNAReadClusterSet::new(read_names_map.clone(), chromosome_names_map.clone());
        cluster_set.add_cluster(RNAReadCluster::new(
            0,
            read_names_map.right_values().copied().collect(),
            Vec::new(),
            variant_calls,
            HashMap::new(),
            HashSet::new()
        ));
        cluster_set
    };
    let run_with = |cluster_set: RNAReadClusterSet, trim_transcript_ends: bool| -> HashMap<Box<str>, Box<str>> {
        let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_transcript_ends, ..CorrectRNAReadsOptions::default() };
        let directory = tempdir().unwrap();
        let output_fastq_file = directory.path().join("corrected.fastq.gz");
        let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
        correct_rna_reads(&bam_file, output_fastq_file, &cluster_set, Some(&annotator), &options, 1, 1000).unwrap();
        open_fastq_reader(output_fastq_file)
            .records()
            .map(|record| {
                let record = record.unwrap();
                (
                    std::str::from_utf8(record.name()).unwrap().into(),
                    std::str::from_utf8(record.sequence()).unwrap().into()
                )
            })
            .collect()
    };
    let run = |cluster_set: RNAReadClusterSet| run_with(cluster_set, true);

    // The reads as sequenced, from the generator's own pieces.
    let tail: String = format!("{}{}{}", &EXON_A[20..60], EXON_B, &INTRON_2[..15]);
    let head: String = format!("{}{}{}", &INTRON_1[29..40], EXON_B, EXON_C);
    let last: String = format!("{}{}{}", EXON_B, EXON_C, &DOWNSTREAM[..20]);

    // Each overhang past an internal exon goes, flush to the exon boundary: the tail loses its
    // last 15 bases and the head its first 11. Exon C is the transcript's last exon, so that
    // end is kept.
    let trimmed = run(build(Vec::new()));
    assert_eq!(trimmed.len(), 3);
    assert_eq!(&*trimmed["tail_overhang"], &tail[..70]);
    assert_eq!(&*trimmed["tail_overhang"], format!("{}{}", &EXON_A[20..60], EXON_B));
    assert_eq!(&*trimmed["head_overhang"], &head[11..]);
    assert_eq!(&*trimmed["head_overhang"], format!("{}{}", EXON_B, EXON_C));
    assert_eq!(&*trimmed["last_exon_overhang"], last);

    // Trimming is off by default: such an end is also what the clusterer counts as a retained
    // intron, so every read comes back as sequenced.
    let untrimmed = run_with(build(Vec::new()), false);
    assert_eq!(&*untrimmed["tail_overhang"], tail);
    assert_eq!(&*untrimmed["head_overhang"], head);
    assert_eq!(&*untrimmed["last_exon_overhang"], last);

    // A called SNV at 140 lies inside the tail's removed span (131..145): that end is kept,
    // while the head, whose span it does not touch, is still cut. The called allele is then
    // applied: the read's own base leaves the tail as sequenced, and another allele rewrites
    // read position 79, which the tail places at 140.
    let snv = |allele: &str| -> GraphOperation {
        GraphOperation::new(
            0, 139, Strand::Forward, GraphOperationType::Downstream,
            0, 141, Strand::Forward, GraphOperationType::Upstream,
            allele.into(), VariantType::SingleNucleotideVariant
        )
    };
    let licensed = run(build(vec![snv(&INTRON_2[9..10])]));
    assert_eq!(&*licensed["tail_overhang"], tail);
    assert_eq!(licensed["head_overhang"], trimmed["head_overhang"]);

    let rewritten = run(build(vec![snv("T")]));
    let mut expected: Vec<u8> = tail.clone().into_bytes();
    assert_ne!(expected[79], b'T');
    expected[79] = b'T';
    assert_eq!(rewritten["tail_overhang"].as_bytes(), expected);
    assert_eq!(rewritten["head_overhang"], trimmed["head_overhang"]);
}
