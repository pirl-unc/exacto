use bimap::BiMap;
use exacto_core::prelude::*;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_reference_transcript_matching_identify_reference_transcript_matches_1() {
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let mut chromosome_names_map: BiMap<Box<str>,u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);

    // CRK - ENST00000574295.1
    let exon_1: TranscriptModelExon = TranscriptModelExon::new(
        0,
        1455877,
        1456149,
        Strand::Reverse,
        1,
        0,
        271
    );
    let exon_2: TranscriptModelExon = TranscriptModelExon::new(
        0,
        1436999,
        1437155,
        Strand::Reverse,
        2,
        272,
        428
    );
    let exon_3: TranscriptModelExon = TranscriptModelExon::new(
        0,
        1422226,
        1423185,
        Strand::Reverse,
        3,
        429,
        1388
    );
    let mut exons: Vec<TranscriptModelExon> = Vec::new();
    exons.push(exon_1);
    exons.push(exon_2);
    exons.push(exon_3);

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &exons,
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    assert_eq!(reference_transcript_matches.len(), 3);
    assert!(reference_transcript_matches.get(0).unwrap().score > 0.9f32);
    assert_eq!(reference_transcript_matches.get(0).unwrap().reference_transcript_id, "ENST00000574295.1".into());
    assert_eq!(reference_transcript_matches.get(1).unwrap().reference_transcript_id, "ENST00000398970.5".into());
    assert_eq!(reference_transcript_matches.get(2).unwrap().reference_transcript_id, "ENST00000300574.3".into());
}
