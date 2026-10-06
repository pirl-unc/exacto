use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_gencode_1() {
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-core/gencode.v41.annotation.tp53.gtf");
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

#[test]
fn test_gencode_ensembl_biotypes_and_unknown_strand() {
    // Ensembl names the types gene_biotype and transcript_biotype, and a GTF writes "." for a
    // feature without a strand (StringTie, for single-exon transcripts).
    let temp_dir = tempfile::TempDir::new().unwrap();
    let gtf_file: String = temp_dir.path().join("ensembl.gtf").to_str().unwrap().to_string();
    let gtf: &str = "\
chrS\ttest\tgene\t1\t1000\t.\t.\t.\tgene_id \"G1\"; gene_name \"A\"; gene_biotype \"protein_coding\"; level 2;
chrS\ttest\ttranscript\t1\t1000\t.\t.\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_biotype \"protein_coding\"; level 2;
chrS\ttest\texon\t1\t1000\t.\t.\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 1; exon_id \"E1\"; level 2;
";
    fs::write(&gtf_file, gtf).unwrap();
    let gencode: Gencode = Gencode::new(
        &gtf_file,
        "hg38",
        "test",
        Some(HashSet::from(["protein_coding"])),
        None,
        Some(HashSet::from(["protein_coding"])),
        None
    );
    let gene: &Gene = gencode.get_gene("G1").unwrap();
    assert_eq!(gene.gene_type, "protein_coding".into());
    assert_eq!(gene.strand, Strand::Unknown);
    let transcript: &Transcript = gencode.get_transcript("T1").unwrap();
    assert_eq!(transcript.transcript_type, "protein_coding".into());
    assert_eq!(transcript.exons.len(), 1);
}

#[test]
#[should_panic(expected = "Loaded 0 genes, 0 transcripts and 0 exons")]
fn test_gencode_no_gene_rows() {
    // Without gene rows no transcript can be placed; the loader stops instead of returning an
    // empty annotation.
    let temp_dir = tempfile::TempDir::new().unwrap();
    let gtf_file: String = temp_dir.path().join("no_genes.gtf").to_str().unwrap().to_string();
    let gtf: &str = "\
chrS\ttest\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_type \"protein_coding\"; level 2;
chrS\ttest\texon\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 1; exon_id \"E1\"; level 2;
";
    fs::write(&gtf_file, gtf).unwrap();
    Gencode::new_with_defaults(&gtf_file, "hg38", "test");
}

#[test]
#[should_panic(expected = "Exon without exon_id")]
fn test_gencode_exon_without_exon_id() {
    // RefSeq and StringTie write no exon_id: every exon of a transcript would share one key and
    // only the last would be kept, so the loader stops.
    let temp_dir = tempfile::TempDir::new().unwrap();
    let gtf_file: String = temp_dir.path().join("no_exon_id.gtf").to_str().unwrap().to_string();
    let gtf: &str = "\
chrS\ttest\tgene\t1\t500\t.\t+\t.\tgene_id \"G1\"; gene_type \"protein_coding\"; level 2;
chrS\ttest\ttranscript\t1\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_type \"protein_coding\"; level 2;
chrS\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 1; level 2;
chrS\ttest\texon\t201\t300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 2; level 2;
chrS\ttest\texon\t401\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number 3; level 2;
";
    fs::write(&gtf_file, gtf).unwrap();
    Gencode::new_with_defaults(&gtf_file, "hg38", "test");
}
