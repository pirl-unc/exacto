use std::collections::HashSet;

use crate::prelude::*;


#[test]
fn test_gene_1() {
    let mut gene: Gene = Gene::new(
        "ENSG001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        "gene",
        1,
        "protein_coding"
    );

    let transcript: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

    gene.add_transcript(transcript);

    assert_eq!(&*gene.get_transcript("ENST001").unwrap().transcript_id, "ENST001");
    assert_eq!(gene.get_transcript_ids().len(), 1);
}

#[test]
fn test_gene_2() {
    let mut gene_1: Gene = Gene::new(
        "ENSG001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        "gene",
        1,
        "protein_coding"
    );

    let transcript: Transcript = Transcript::new(
        "ENSG001",
        "ENST001",
        "havana",
        "chr1",
        1001,
        2000,
        Strand::Forward,
        1,
        "transcript",
        "protein_coding",
        "1",
        HashSet::new()
    );

    gene_1.add_transcript(transcript);

    let gene_2: Gene = gene_1.clone();
    
    assert_eq!(gene_1, gene_2);
}
