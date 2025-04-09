use std::fs;
use std::path::Path;

use crate::structs::gene_annotation::gene_annotator::GeneAnnotator;
use crate::structs::gene_annotation::gene::Gene;
use crate::structs::gene_annotation::transcript::Transcript;
use crate::structs::gene_annotation::gencode::Gencode;
use crate::structs::gene_annotation::intron::Intron;


#[test]
fn test_gencode_1() {
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let mut gencode: Gencode = Gencode::new(gtf_file, "hg38");

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
fn test_get_introns_1() {
    let gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gtf_file: &str = gtf_full_path.to_str().unwrap();
    let mut gencode: Gencode = Gencode::new(gtf_file, "hg38");
    let transcript: &Transcript = gencode.get_transcript("ENST00000269305.9").unwrap();
    let introns: Vec<Intron> = transcript.get_introns();
    assert!(introns.len() == 10);
    assert!(introns[0].start == 7676623);
    assert!(introns[0].end == 7687376);
}
