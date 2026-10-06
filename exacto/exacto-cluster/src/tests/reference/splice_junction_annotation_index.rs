use exacto_core::prelude::Gencode;
use std::fs;
use std::path::Path;

use super::*;


/// The index of the annotation of chr17 and chr18 (GENCODE v41): 19,211 reference transcripts,
/// 17,349 of them with an intron, 101,199 introns over 29,588 junctions and 29,220 pairs of
/// consecutive junctions.
///
/// The first intron of TP53 (chr17, reverse strand), 7,676,623-7,687,376, belongs to one gene
/// and to 11 reference transcripts. The second intron, 7,676,404-7,676,520, follows it in
/// transcript order and does not precede it.
#[test]
fn splice_junction_annotation_index_returns_matches() {
    let gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([
        (Box::<str>::from("chr17"), 0u16),
        (Box::<str>::from("chr18"), 1u16)
    ]);

    let index: SpliceJunctionAnnotationIndex = SpliceJunctionAnnotationIndex::new(&gene_annotator, &chromosome_names_map);

    assert_eq!(index.transcripts.len(), 17_349);
    assert_eq!(index.junctions_to_genes.len(), 29_588);
    assert_eq!(index.junctions_to_transcripts.len(), 29_588);
    assert_eq!(index.junctions_to_transcripts.values().map(|transcripts| transcripts.len()).sum::<usize>(), 101_199);
    assert_eq!(index.annotated_transitions.len(), 29_220);
    assert_eq!(index.gene_spans[&0].len(), 3_103);
    assert_eq!(index.gene_spans[&1].len(), 1_260);

    assert_eq!(index.gene_spans[&0]["ENSG00000141510.18"], (7_661_779, 7_687_538));
    assert_eq!(
        index.junctions_to_genes[&(0, 7_676_623, 7_687_376)],
        BTreeSet::from([Box::<str>::from("ENSG00000141510.18")])
    );
    let mut reference_gene_transcript_ids: Vec<(&str, &str)> = index.junctions_to_transcripts[&(0, 7_676_623, 7_687_376, Strand::Reverse)]
        .iter()
        .map(|&reference| (&*index.transcripts[reference as usize].0, &*index.transcripts[reference as usize].1))
        .collect();
    reference_gene_transcript_ids.sort();
    assert_eq!(
        reference_gene_transcript_ids,
        vec![
            ("ENSG00000141510.18", "ENST00000269305.9"),
            ("ENSG00000141510.18", "ENST00000420246.6"),
            ("ENSG00000141510.18", "ENST00000455263.6"),
            ("ENSG00000141510.18", "ENST00000505014.5"),
            ("ENSG00000141510.18", "ENST00000514944.5"),
            ("ENSG00000141510.18", "ENST00000604348.5"),
            ("ENSG00000141510.18", "ENST00000610292.4"),
            ("ENSG00000141510.18", "ENST00000610538.4"),
            ("ENSG00000141510.18", "ENST00000620739.4"),
            ("ENSG00000141510.18", "ENST00000622645.4"),
            ("ENSG00000141510.18", "ENST00000635293.1")
        ]
    );

    let first_intron: SpliceJunction = SpliceJunction::new(0, 0, 7_687_376, 7_676_623, Strand::Reverse, Strand::Reverse);
    let second_intron: SpliceJunction = SpliceJunction::new(0, 0, 7_676_520, 7_676_404, Strand::Reverse, Strand::Reverse);
    assert!(index.annotated_transitions.contains(&(first_intron.clone(), second_intron.clone())));
    assert!(!index.annotated_transitions.contains(&(second_intron, first_intron)));
}
