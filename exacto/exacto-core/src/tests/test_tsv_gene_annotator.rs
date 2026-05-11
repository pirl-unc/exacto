use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_tsv_gene_annotator_1() {
    let tsv_path = Path::new("src/tests/data/tsv/sample_gene_annotations.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();
    let tsv_gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(tsv_file, "custom", "v1");

    assert!(tsv_gene_annotator.get_assembly() == "custom");
    assert!(tsv_gene_annotator.get_version() == "v1");
    assert!(tsv_gene_annotator.get_gene_ids_at_locus("chrA", 25).len() == 1);
    assert!(tsv_gene_annotator.get_transcripts().len() == 4);
    assert!(tsv_gene_annotator.get_transcript_ids_overlapping_region("chrA", 1, 5).len() == 2);
    assert!(tsv_gene_annotator.get_transcript_ids_overlapping_region("chrA", 46, 50).len() == 1);

    let gene: &Gene = tsv_gene_annotator.get_gene("g1").unwrap();
    assert!(gene.gene_id == "g1".into());
    assert!(gene.start == 1);
    assert!(gene.end == 50);

    let transcript: &Transcript = tsv_gene_annotator.get_transcript("t1").unwrap();
    assert!(transcript.transcript_id == "t1".into());
    assert!(transcript.start == 1);
    assert!(transcript.end == 25);

    let transcript: &Transcript = tsv_gene_annotator.get_transcript("t2").unwrap();
    assert!(transcript.transcript_id == "t2".into());
    assert!(transcript.start == 1);
    assert!(transcript.end == 50);
}

#[test]
fn test_tsv_gene_annotator_2() {
    let tsv_path = Path::new("src/tests/data/tsv/sample_gene_annotations.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();
    let tsv_gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(tsv_file, "custom", "v1");

    let transcript: &Transcript = tsv_gene_annotator.get_transcript("t1").unwrap();
    let introns: Vec<Intron> = transcript.get_introns();

    assert!(introns.len() == 2);
    assert!(introns[0].start == 6);
    assert!(introns[0].end == 10);
    assert!(introns[1].start == 16);
    assert!(introns[1].end == 20);

    assert_eq!(tsv_gene_annotator.get_exon_ids_overlapping_region("chrA",21, 25).len(), 2);
    assert_eq!(tsv_gene_annotator.get_exon("t1","e1").unwrap().exon_number, 1);
    assert_eq!(tsv_gene_annotator.get_exons().len(), 13);
}
