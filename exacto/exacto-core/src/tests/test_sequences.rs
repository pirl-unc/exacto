use std::collections::{HashMap, HashSet};
use crate::prelude::*;


#[test]
fn test_find_kmers_1() {
    let value: String = "AAATTTCCCAAA".to_string();
    let kmers: HashMap<Box<str>,Vec<usize>> = find_kmers(value.as_str(), 3);
    match kmers.get("AAA") {
        Some(values) => {
            assert_eq!(values[0], 0);
            assert_eq!(values[1], 9);
        }
        None => {
            panic!("Unexpected error.")
        }
    }
}

#[test]
fn test_is_valid_nucleotide_sequence_1() {
    assert!(is_valid_nucleotide_sequence("ATCGACGactg"));
    assert!(is_valid_nucleotide_sequence("AUCGACGaucg"));
    assert!(is_valid_nucleotide_sequence("AUCGTCG") == false);
}

#[test]
fn test_reverse_complement_1() {
    let seq: String = "ATCG".to_string();
    let reverse_complement: Box<str> = reverse_complement(seq.as_str());
    assert_eq!(reverse_complement, "CGAT".into());
}

#[test]
fn test_reverse_complement_2() {
    let seq: String = "cgat".to_string();
    let reverse_complement: Box<str> = reverse_complement(seq.as_str());
    assert_eq!(reverse_complement, "atcg".into());
}

#[test]
fn test_translate_1() {
    let rna_sequence: String = "AUGUAG".to_string();
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    assert_eq!(peptides.len(), 1);
    assert_eq!(peptides[0].0, "M*".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 5);
    assert_eq!(peptides[0].3, 2);
}

#[test]
fn test_translate_2() {
    let rna_sequence: String = "AUGAGUAUCAUCAACUUUGAAAAACUCUAG".to_string();
    let peptides: Vec<(Box<str>, usize, usize, usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides[0].0, "MSIINFEKL*".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 29);
    assert_eq!(peptides[0].3, 10);
    assert_eq!(peptides[1].0, "LKNS".into());
    assert_eq!(peptides[1].1, 16);
    assert_eq!(peptides[1].2, 27);
    assert_eq!(peptides[1].3, 4);
}

#[test]
fn test_translate_3() {
    let rna_sequence: String = "AAUGAGUAUCAUCAACUUUGAAAAACUCUAGAAAAAAAUGUGUUGUUGUAUCAUCAACUUUGAAAAACUCUAG".to_string();
    let peptides: Vec<(Box<str>, usize, usize, usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }
    assert_eq!(peptides.len(), 7);
    assert_eq!(peptides_set.contains("LKNS"), true);
    assert_eq!(peptides_set.contains("LYHQL*"), true);
    assert_eq!(peptides_set.contains("LLYHQL*"), true);
    assert_eq!(peptides_set.contains("VLLYHQL*"), true);
    assert_eq!(peptides_set.contains("MCCCIINFEKL*"), true);
    assert_eq!(peptides_set.contains("LKNSRKKCVVVSSTLKNS"), true);
    assert_eq!(peptides_set.contains("MSIINFEKL*"), true);
}

#[test]
fn test_translate_4() {
    let rna_sequence: String = "AUGAUUUGCCAUAUCGGGGCGAAC".to_string();
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides[0].0, "MICHIGAN".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 23);
    assert_eq!(peptides[0].3, 8);
    assert_eq!(peptides[1].0, "LPYRGE".into());
    assert_eq!(peptides[1].1, 5);
    assert_eq!(peptides[1].2, 22);
    assert_eq!(peptides[1].3, 6);
}

#[test]
fn test_translate_5() {
    let rna_sequence: String = "AUGAUUUGCCAUAUCGGGGCGAACUGAAUGAUUUGCCAUAUCGGGGCGAAC".to_string();
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }
    assert_eq!(peptides.len(), 5);
    assert_eq!(peptides_set.contains("LNDLPYRGE"), true);
    assert_eq!(peptides_set.contains("LPYRGELNDLPYRGE"), true);
    assert_eq!(peptides_set.contains("MICHIGAN*"), true);
    assert_eq!(peptides_set.contains("LPYRGE"), true);
    assert_eq!(peptides_set.contains("MICHIGAN"), true);
}

#[test]
fn test_translate_6() {
    let rna_sequence: String = "AUUAUUAUUAUUAUUAUUAUUAUUAUU".to_string();
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    assert!(peptides.is_empty());
}

#[test]
fn test_translate_7() {
    let rna_sequence: String = "ATGGGGCCCATGCCTTAG".to_string();
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }
    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides_set.contains("MGPMP*"), true);
    assert_eq!(peptides_set.contains("MP*"), true);
}
