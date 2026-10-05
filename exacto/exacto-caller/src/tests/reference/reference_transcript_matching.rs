use bimap::BiMap;
use std::fs;
use std::path::Path;
use exacto_core::prelude::Gencode;

use crate::prelude::*;

use super::*;


fn get_gencode() -> Gencode {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    Gencode::new_with_defaults(
        gencode_gtf_full_path.to_str().unwrap(),
        "hg38",
        "v41",
    )
}

fn get_chromosome_names_map() -> BiMap<Box<str>, u16> {
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    chromosome_names_map
}


#[test]
fn reference_transcript_matching_returns_correct_transcript_match() {
    let gene_annotator: Gencode = get_gencode();
    let chromosome_names_map: BiMap<Box<str>, u16> = get_chromosome_names_map();

    // CRK gene
    let exons: Vec<TranscriptModelExon> = vec![
        TranscriptModelExon::new(0, 1455877, 1456149, Strand::Reverse, 1, 0, 271),
        TranscriptModelExon::new(0, 1436999, 1437155, Strand::Reverse, 2, 272, 428),
        TranscriptModelExon::new(0, 1422226, 1423185, Strand::Reverse, 3, 429, 1388),
    ];
    let introns: Vec<TranscriptModelSpliceJunction> = vec![
        // exon 1 (1455877..) <-> exon 2 (..1437155): 1437156 ..= 1455876
        TranscriptModelSpliceJunction::new(
            0,
            0,
            1437156,
            1455876,
            Strand::Reverse,
            Strand::Reverse,
            1,
            271,
            272),
        // exon 2 (1436999..) <-> exon 3 (..1423185): 1423186 ..= 1436998
        TranscriptModelSpliceJunction::new(
            0,
            0,
            1423186,
            1436998,
            Strand::Reverse,
            Strand::Reverse,
            2,
            428,
            429),
    ];

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> =
        identify_reference_transcript_matches(
            &exons,
            &introns,
            &gene_annotator,
            &chromosome_names_map,
        );

    assert_eq!(reference_transcript_matches.len(), 1);

    let best_match: &ReferenceTranscriptMatch = &reference_transcript_matches[0];
    assert_eq!(best_match.get_reference_gene_id(), "ENSG00000167193.8");
    assert_eq!(best_match.get_reference_transcript_id(), "ENST00000574295.1");

    assert_eq!(best_match.num_splice_junction_matches(), 1);
    assert_eq!(best_match.num_query_only_bases(), 0);
    assert_eq!(best_match.num_reference_only_bases(), 1);
    assert_eq!(best_match.num_overlapping_bases(), 1390);
}


#[test]
fn reference_transcript_matching_returns_protein_coding_over_short_nmd_transcript() {
    let gene_annotator: Gencode = get_gencode();
    let chromosome_names_map: BiMap<Box<str>, u16> = get_chromosome_names_map();

    // SERPINF1 gene
    let exons: Vec<TranscriptModelExon> = vec![
        TranscriptModelExon::new(0, 1762104, 1762113, Strand::Forward, 1, 0, 9),
        TranscriptModelExon::new(0, 1766903, 1766994, Strand::Forward, 2, 10, 101),
        TranscriptModelExon::new(0, 1769852, 1769861, Strand::Forward, 3, 102, 111),
    ];
    let introns: Vec<TranscriptModelSpliceJunction> = vec![
        TranscriptModelSpliceJunction::new(
            0,
            0,
            1762114,
            1766902,
            Strand::Forward,
            Strand::Forward,
            1,
            9,
            10),
        TranscriptModelSpliceJunction::new(
            0,
            0,
            1766995,
            1769851,
            Strand::Forward,
            Strand::Forward,
            2,
            101,
            102),
    ];

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> =
        identify_reference_transcript_matches(
            &exons,
            &introns,
            &gene_annotator,
            &chromosome_names_map,
        );

    let serpinf1_match: &ReferenceTranscriptMatch = reference_transcript_matches
        .iter()
        .find(|m| m.get_reference_gene_id() == "ENSG00000132386.11")
        .expect("expected a SERPINF1 match");

    assert_eq!(serpinf1_match.num_splice_junction_matches(), 2);
    assert_eq!(
        serpinf1_match.get_reference_transcript_id(),
        "ENST00000254722.9",
        "the protein_coding isoform must beat the shorter NMD isoform ENST00000573770.5 \
         despite its larger reference-only base count"
    );
}



/// The matches of a query that overlaps two genes come best first, whatever order the genes
/// are visited in.
///
/// GLTPD2 and PSMB6 lie on the forward strand of chr17, 5.6 kb apart. Both queries are
/// read-throughs: exons of GLTPD2, a junction into PSMB6, exons of PSMB6.
///
///   Query 1   GLTPD2 exons 1-4 (3 of its junctions), PSMB6 exons 2-6 (4 of its junctions)
///   Query 2   GLTPD2 exons 1-4 (3 of its junctions), PSMB6 exons 2-3 (1 of its junctions)
///
/// The match with more splice junctions is the better one: PSMB6 in query 1, GLTPD2 in
/// query 2. The first exon of GLTPD2 lies inside ENSG00000280254, a gene of one exon, whose
/// transcript matches no junction and comes last in both.
#[test]
fn identify_reference_transcript_matches_returns_matches_best_first() {
    let gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([
        (Box::<str>::from("chr17"), 0u16)
    ]);

    // (exons of the query, (reference transcript ID, number of splice junction matches))
    for (exon_spans, matches) in [
        (
            vec![
                (4_788_964u32, 4_789_117u32), (4_789_226, 4_789_290), (4_789_511, 4_789_677), (4_789_759, 4_790_000),
                (4_796_728, 4_796_795), (4_797_438, 4_797_569), (4_797_682, 4_797_811), (4_798_009, 4_798_153), (4_798_280, 4_798_495)
            ],
            vec![("ENST00000270586.8", 4usize), ("ENST00000331264.8", 3usize), ("ENST00000623798.1", 0usize)]
        ),
        (
            vec![
                (4_788_964, 4_789_117), (4_789_226, 4_789_290), (4_789_511, 4_789_677), (4_789_759, 4_790_000),
                (4_796_728, 4_796_795), (4_797_438, 4_797_569)
            ],
            vec![("ENST00000331264.8", 3), ("ENST00000270586.8", 1), ("ENST00000623798.1", 0)]
        )
    ] {
        let mut exons: Vec<TranscriptModelExon> = Vec::new();
        let mut splice_junctions: Vec<TranscriptModelSpliceJunction> = Vec::new();
        let mut read_position: u32 = 0;
        for (i, &(exon_start, exon_end)) in exon_spans.iter().enumerate() {
            if i > 0 {
                splice_junctions.push(TranscriptModelSpliceJunction::new(
                    0,
                    0,
                    exon_spans[i - 1].1 + 1,
                    exon_start - 1,
                    Strand::Forward,
                    Strand::Forward,
                    i as u16,
                    read_position - 1,
                    read_position
                ));
            }
            exons.push(TranscriptModelExon::new(
                0,
                exon_start,
                exon_end,
                Strand::Forward,
                (i + 1) as u16,
                read_position,
                read_position + (exon_end - exon_start)
            ));
            read_position += exon_end - exon_start + 1;
        }

        // The genes are visited in the order of a hash set: every call may visit them in
        // another order.
        for _ in 0..20 {
            let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
                &exons,
                &splice_junctions,
                &gene_annotator,
                &chromosome_names_map
            );

            assert_eq!(
                reference_transcript_matches
                    .iter()
                    .map(|reference_transcript_match| (
                        reference_transcript_match.get_reference_transcript_id(),
                        reference_transcript_match.num_splice_junction_matches()
                    ))
                    .collect::<Vec<(&str, usize)>>(),
                matches
            );
        }
    }
}

#[test]
fn identify_reference_transcript_matches_gives_the_gene_name_and_the_gene_id() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gencode_gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);

    // Three exons of CRK (ENSG00000167193.8).
    let exons: Vec<TranscriptModelExon> = vec![
        TranscriptModelExon::new(0, 1455877, 1456149, Strand::Reverse, 1, 0, 271),
        TranscriptModelExon::new(0, 1436999, 1437155, Strand::Reverse, 2, 272, 428),
        TranscriptModelExon::new(0, 1422226, 1423185, Strand::Reverse, 3, 429, 1388),
    ];
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = vec![
        TranscriptModelSpliceJunction::new(0, 0, 1437156, 1455876, Strand::Reverse, Strand::Reverse, 1, 271, 272),
        TranscriptModelSpliceJunction::new(0, 0, 1423186, 1436998, Strand::Reverse, Strand::Reverse, 2, 428, 429),
    ];

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &exons,
        &splice_junctions,
        &gene_annotator,
        &chromosome_names_map
    );

    assert_eq!(reference_transcript_matches[0].get_reference_gene_id(), "ENSG00000167193.8");
    assert_eq!(reference_transcript_matches[0].get_reference_gene_name(), "CRK");
}
