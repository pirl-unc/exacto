use std::collections::HashMap;

use crate::common::constants::*;
use crate::common::sequences::*;


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
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
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
    let peptides: Vec<(Box<str>,usize,usize,usize)> = translate(
        rna_sequence.as_str(),
        START_CODONS.iter().map(|c| &**c).collect()
    );
    assert_eq!(peptides.len(), 4);
    assert_eq!(peptides[0].0, "VLLYHQL*".into());
    assert_eq!(peptides[0].1, 39);
    assert_eq!(peptides[0].2, 62);
    assert_eq!(peptides[0].3, 8);
    assert_eq!(peptides[1].0, "MSIINFEKL*".into());
    assert_eq!(peptides[1].1, 1);
    assert_eq!(peptides[1].2, 30);
    assert_eq!(peptides[1].3, 10);
    assert_eq!(peptides[2].0, "MCCCIINFEKL*".into());
    assert_eq!(peptides[2].1, 37);
    assert_eq!(peptides[2].2, 72);
    assert_eq!(peptides[2].3, 12);
    assert_eq!(peptides[3].0, "LKNSRKKCVVVSSTLKNS".into());
    assert_eq!(peptides[3].1, 17);
    assert_eq!(peptides[3].2, 70);
    assert_eq!(peptides[3].3, 18);
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
    assert_eq!(peptides.len(), 3);
    assert_eq!(peptides[0].0, "MICHIGAN*".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 26);
    assert_eq!(peptides[0].3, 9);
    assert_eq!(peptides[1].0, "MICHIGAN".into());
    assert_eq!(peptides[1].1, 27);
    assert_eq!(peptides[1].2, 50);
    assert_eq!(peptides[1].3, 8);
    assert_eq!(peptides[2].0, "LPYRGELNDLPYRGE".into());
    assert_eq!(peptides[2].1, 5);
    assert_eq!(peptides[2].2, 49);
    assert_eq!(peptides[2].3, 15);
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

