use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_gencode_1() {
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let mut gencode: Gencode = Gencode::new_with_defaults(gtf_file,"hg38", "v41");

    let gene_ids: Vec<Box<str>> = gencode.get_gene_ids_overlapping_region("chr17", 7600000, 7700000);
    assert!(gene_ids.len() == 6);

    let gene: &Gene = gencode.get_gene("ENSG00000141510.18").unwrap();
    assert!(gene.gene_id == "ENSG00000141510.18".into());
    assert!(gene.start == 7661779);
    assert!(gene.end == 7687538);
    assert!(gene.gene_type == "protein_coding".into());

    let transcript: &Transcript = gencode.get_transcript("ENST00000413465.6").unwrap();
    assert!(transcript.transcript_id == "ENST00000413465.6".into());
    assert!(transcript.start == 7661779);
    assert!(transcript.end == 7676594);
    assert!(transcript.transcript_type == "protein_coding".into());
}

#[test]
fn test_gencode_2() {
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let mut gencode: Gencode = Gencode::new_with_defaults(gtf_file,"hg38", "v41");
    let transcript: &Transcript = gencode.get_transcript("ENST00000269305.9").unwrap();
    let introns: Vec<Intron> = transcript.get_introns();
    assert!(introns.len() == 10);
    assert!(introns[0].start == 7676623);
    assert!(introns[0].end == 7687376);
}

#[test]
fn test_gencode_3() {
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.tp53.gtf");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let mut gencode: Gencode = Gencode::new_with_defaults(gtf_file,"hg38", "v41");
    let transcript: &Transcript = gencode.get_transcript("ENST00000269305.9").unwrap();
    let introns: Vec<Intron> = transcript.get_introns();
    assert_eq!(gencode.get_assembly(), "hg38");
    assert_eq!(gencode.get_version(), "v41");
    assert_eq!(gencode.get_gene_ids_at_locus("chr17", 7676623).len(), 1);
    assert_eq!(gencode.get_transcript_ids_overlapping_region("chr17",7_668_000, 7_689_000).len(), 27);
    assert_eq!(gencode.get_transcript_ids_overlapping_region("chr18",7_668_000, 7_689_000).len(), 0);
    assert_eq!(gencode.get_exon_ids_overlapping_region("chr17",7_670_600, 7_670_700).len(), 16);
    assert_eq!(gencode.get_exon("ENST00000503591.1","ENSE00002419584.1").unwrap().exon_number, 4);
    assert_eq!(gencode.get_genes().len(), 1);
    assert_eq!(gencode.get_transcripts().len(), 27);
    assert_eq!(gencode.get_exons().len(), 214);
    assert_eq!(introns.len(), 10);
    assert_eq!(introns[0].start, 7676623);
    assert_eq!(introns[0].end, 7687376);
}
