use std::collections::{HashMap, HashSet};

use crate::prelude::*;


#[test]
fn test_find_kmers_1() {
    let value: String = "AAATTTCCCAAA".to_string();
    let kmers: HashMap<Box<str>, Vec<u32>> = find_kmers(value.as_str(), 3);
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
fn test_find_substring_positions_1() {
    let value: String = "helloworldhelloworld".to_string();
    let positions: Vec<u32> = find_substring_positions(value.as_str(), "hello");
    assert_eq!(positions, vec![0, 10]);
}

#[test]
fn test_get_dinucleotide_repeat_motif_1() {
    assert_ne!(get_dinucleotide_repeat_motif("AC"), None);
    assert!(get_dinucleotide_repeat_motif("AC").unwrap() == "AC".to_string());
    assert!(get_dinucleotide_repeat_motif("CAC").unwrap() == "CA".to_string());
}

#[test]
fn test_is_homopolymer_sequence_1() {
    assert!(is_homopolymer_sequence("A") == true);
    assert!(is_homopolymer_sequence("AAAAAA") == true);
    assert!(is_homopolymer_sequence("AAAaaa") == true);
    assert!(is_homopolymer_sequence("cCcccC") == true);
    assert!(is_homopolymer_sequence("Gggggg") == true);
    assert!(is_homopolymer_sequence("AAATtt") == false);
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
fn test_reverse_complement_iupac() {
    // Every ambiguity code a BAM `SEQ` can hold, so a reverse-strand read carrying one
    // still turns back into read orientation.
    assert_eq!(reverse_complement("ACGTRYKMSWBDHVN"), "NBDHVWSKMRYACGT".into());
    assert_eq!(reverse_complement("rykmswbdhvn"), "nbdhvwskmry".into());
}

#[test]
fn test_reverse_string_1() {
    assert!(reverse_string("abcd") == "dcba");
}


#[test]
fn test_translate_1() {
    let rna_sequence: &str = "AUGUAG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);

    assert_eq!(peptides.len(), 1);
    assert_eq!(peptides[0].0, "M*".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 5);
    assert_eq!(peptides[0].3, 2);
}

#[test]
fn test_translate_2() {
    let rna_sequence: &str = "AUGAGUAUCAUCAACUUUGAAAAACUCUAG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);

    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides[0].0, "MSIINFEKL*".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 29);
    assert_eq!(peptides[0].3, 10);
    assert_eq!(peptides[1].0, "MKNS".into());
    assert_eq!(peptides[1].1, 16);
    assert_eq!(peptides[1].2, 27);
    assert_eq!(peptides[1].3, 4);
}

#[test]
fn test_translate_3() {
    let rna_sequence: &str = "AAUGAGUAUCAUCAACUUUGAAAAACUCUAGAAAAAAAUGUGUUGUUGUAUCAUCAACUUUGAAAAACUCUAG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }

    assert_eq!(peptides.len(), 7);
    assert_eq!(peptides_set.contains("MKNS"), true);
    assert_eq!(peptides_set.contains("MYHQL*"), true);
    assert_eq!(peptides_set.contains("MLYHQL*"), true);
    assert_eq!(peptides_set.contains("MLLYHQL*"), true);
    assert_eq!(peptides_set.contains("MCCCIINFEKL*"), true);
    assert_eq!(peptides_set.contains("MKNSRKKCVVVSSTLKNS"), true);
    assert_eq!(peptides_set.contains("MSIINFEKL*"), true);
}

#[test]
fn test_translate_4() {
    let rna_sequence: &str = "AUGAUUUGCCAUAUCGGGGCGAAC";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);

    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides[0].0, "MICHIGAN".into());
    assert_eq!(peptides[0].1, 0);
    assert_eq!(peptides[0].2, 23);
    assert_eq!(peptides[0].3, 8);
    assert_eq!(peptides[1].0, "MPYRGE".into());
    assert_eq!(peptides[1].1, 5);
    assert_eq!(peptides[1].2, 22);
    assert_eq!(peptides[1].3, 6);
}

#[test]
fn test_translate_5() {
    let rna_sequence: &str = "AUGAUUUGCCAUAUCGGGGCGAACUGAAUGAUUUGCCAUAUCGGGGCGAAC";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }

    assert_eq!(peptides.len(), 5);
    assert_eq!(peptides_set.contains("MNDLPYRGE"), true);
    assert_eq!(peptides_set.contains("MPYRGELNDLPYRGE"), true);
    assert_eq!(peptides_set.contains("MICHIGAN*"), true);
    assert_eq!(peptides_set.contains("MPYRGE"), true);
    assert_eq!(peptides_set.contains("MICHIGAN"), true);
}

#[test]
fn test_translate_6() {
    let rna_sequence: &str = "AUUAUUAUUAUUAUUAUUAUUAUUAUU";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);

    assert!(peptides.is_empty());
}

#[test]
fn test_translate_7() {
    let rna_sequence: &str = "ATGGGGCCCATGCCTTAG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let peptides = translate(rna_sequence, &start_codons);
    let mut peptides_set: HashSet<&str> = HashSet::new();
    for peptide in peptides.iter() {
        peptides_set.insert(&peptide.0);
    }

    assert_eq!(peptides.len(), 2);
    assert_eq!(peptides_set.contains("MGPMP*"), true);
    assert_eq!(peptides_set.contains("MP*"), true);
}

#[test]
fn test_translate_8() {
    // An ORF opened at GUG, CUG or UUG begins with methionine, as one opened at AUG does; the
    // same codon inside the ORF keeps its own residue.
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    for (rna_sequence, peptide) in [("GUGAAAGUGUAA", "MKV*"), ("CUGAAACUGUAA", "MKL*"), ("UUGAAAUUGUAA", "MKL*")] {
        let orfs: Vec<(Box<str>, u32, u32, u32)> = identify_open_reading_frames(
            rna_sequence,
            &TranslationStrategy::LongestORF,
            &start_codons
        );
        assert_eq!(orfs, vec![(peptide.into(), 0, 11, 4)]);
    }
}

#[test]
fn test_identify_open_reading_frames_all_orfs_keeps_every_complete_frame() {
    // AUG at 0 closes at 8 (MP*), AUG at 9 closes at 20 (MKK*), AUG at 21 runs off the end.
    let rna_sequence: &str = "AUGCCUUAGAUGAAAAAAUAGAUGGG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let orfs: Vec<(Box<str>, u32, u32, u32)> = identify_open_reading_frames(
        rna_sequence,
        &TranslationStrategy::AllORFs,
        &start_codons
    );
    assert_eq!(orfs, vec![("MP*".into(), 0, 8, 3), ("MKK*".into(), 9, 20, 4)]);
}

#[test]
fn test_identify_open_reading_frames_longest_orf_keeps_one() {
    let rna_sequence: &str = "AUGCCUUAGAUGAAAAAAUAGAUGGG";
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    let orfs: Vec<(Box<str>, u32, u32, u32)> = identify_open_reading_frames(
        rna_sequence,
        &TranslationStrategy::LongestORF,
        &start_codons
    );
    assert_eq!(orfs, vec![("MKK*".into(), 9, 20, 4)]);

    // Equal lengths: the last frame wins, as the translator has always resolved ties.
    let orfs: Vec<(Box<str>, u32, u32, u32)> = identify_open_reading_frames(
        "AUGUAGAUGUAA",
        &TranslationStrategy::LongestORF,
        &start_codons
    );
    assert_eq!(orfs, vec![("M*".into(), 6, 11, 2)]);
}

#[test]
fn test_identify_open_reading_frames_incomplete_frames_are_dropped() {
    let start_codons: HashSet<&str> = START_CODONS.iter().map(|c| c.as_ref()).collect();
    for strategy in [TranslationStrategy::AllORFs, TranslationStrategy::LongestORF] {
        assert!(
            identify_open_reading_frames("AUGGCUGCU", &strategy, &start_codons).is_empty()
        );
    }
}
