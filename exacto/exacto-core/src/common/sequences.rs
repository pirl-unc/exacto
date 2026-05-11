// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use std::collections::{HashMap, HashSet};

use crate::common::constants::CODON_TABLE;


/// Find all k-mers in a sequence.
///
/// Returns:
///
/// * A HashMap where the key is k-mer and the value is a vector of start positions in the sequence.
pub fn find_kmers(sequence: &str, k: u32) -> HashMap<Box<str>, Vec<u32>> {
    let mut kmer_index: HashMap<Box<str>,Vec<u32>> = HashMap::new();
    if sequence.len() < k as usize {
        return kmer_index;
    }
    for i in 0..=sequence.len() - k as usize {
        let kmer: &str = &sequence[i..i + k as usize];
        kmer_index.entry(kmer.to_string().into_boxed_str())
            .or_insert(Vec::new())
            .push(i as u32);
    }
    kmer_index
}

/// Find substring start positions.
pub fn find_substring_positions(text: &str, substring: &str) -> Vec<u32> {
    let mut start_positions: Vec<u32> = Vec::new();
    let substring_length: u32 = substring.len() as u32;
    if substring_length == 0 || substring_length > text.len() as u32 {
        return start_positions;
    }
    for (index, window) in text.as_bytes().windows(substring_length as usize).enumerate() {
        if window == substring.as_bytes() {
            start_positions.push(index as u32);
        }
    }
    start_positions
}

pub fn get_dinucleotide_repeat_motif(seq: &str) -> Option<String> {
    let s = seq.as_bytes();
    let n = s.len();

    if n < 2 {
        return None;
    }

    let a = s[0];
    let b = s[1];

    // Reject homopolymer
    if a == b {
        return None;
    }

    // Check alternating pattern
    for i in 0..n {
        let expected = if (i & 1) == 0 { a } else { b };

        if s[i] != expected {
            return None;
        }
    }

    let motif: String = String::from_utf8(vec![a, b]).ok()?;

    Some(motif)
}

pub fn is_homopolymer_sequence(s: &str) -> bool {
    if s.is_empty() {
        return false;
    }
    let first = s.chars().next().unwrap().to_ascii_uppercase();
    s.chars().all(|c| c.to_ascii_uppercase() == first)
}

/// Check if a sequence is a valid nucleotide (DNA or RNA) sequence.
///
/// Notes:
///
/// * If the sequence is a RNA sequence (i.e. there is a U or u in the sequence), the sequence must not have any T or t.
/// * If the sequence is a DNA sequence (i.e. there is a T or t in the sequence), the sequence must not have any U or u.
pub fn is_valid_nucleotide_sequence(sequence: &str) -> bool {
    if sequence.contains('U') || sequence.contains('u') {
        // RNA sequence: check if all characters are valid RNA nucleotides
        sequence.chars().all(|c| matches!(c, 'A' | 'U' | 'C' | 'G' | 'a' | 'u' | 'c' | 'g'))
    } else if sequence.contains('T') || sequence.contains('t') {
        // DNA sequence: check if all characters are valid DNA nucleotides
        sequence.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G' | 'a' | 't' | 'c' | 'g'))
    } else {
        // Check if all characters are 'A', 'C', 'G' only.
        sequence.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'N' | 'a' | 'c' | 'g' | 'n'))
    }
}

/// Get the reverse complement of a nucleotide sequence.
///
/// Parameters:
///
/// * `sequence` is a nucleotide sequence (DNA or RNA).
///
/// Returns:
///
/// * Reverse complement sequence.
pub fn reverse_complement(sequence: &str) -> Box<str> {
    sequence.chars()
        .rev()
        .map(|nucleotide| match nucleotide {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            'a' => 't',
            't' => 'a',
            'u' => 'a',
            'c' => 'g',
            'g' => 'c',
            'n' => 'n',
            _ => {
                panic!("Invalid nucleotide: {}", nucleotide);
            }
        })
        .collect()
}

pub fn reverse_string(s: &str) -> String {
    s.chars().rev().collect()
}

/// Translate an RNA sequence to all possible ORF peptides.
///
/// Parameters:
///
/// * `rna_sequence` is an RNA sequence consisting of characters \[AUCGaucg\]+.
/// * `start_codons` is a vector of possible start codons (e.g. 'AUG', 'GUG', 'CUG', 'UUG').
///
/// Returns:
///
/// * A vector of tuples (peptide sequence, ORF start, ORF end, peptide length).
pub fn translate(
    rna_sequence: &str,
    start_codons: HashSet<&str>
) -> Vec<(Box<str>, u32, u32, u32)> {
    // Step 1. Convert DNA to RNA
    let rna_sequence_: Box<str> = rna_sequence
        .replace('T', "U")
        .replace('t', "u")
        .into_boxed_str();

    // Step 2. Make sure the start codons are all upper case
    let mut start_codons_: HashSet<Box<str>> = HashSet::new();
    for start_codon in start_codons.iter() {
        start_codons_.insert(start_codon.to_uppercase().into());
    }

    // Step 3. Identify all start codons
    let mut start_codon_positions: Vec<u32> = Vec::new();
    for i in 0..rna_sequence_.len() {
        if i + 2 >= rna_sequence_.len() {
            break;
        }
        let codon: String = rna_sequence_[i..=i+2].to_uppercase();
        if start_codons_.contains(&codon.into_boxed_str()) {
            start_codon_positions.push(i as u32);
        }
    }

    // Step 4. Translate from each start codon
    let mut peptides: Vec<(Box<str>, u32, u32, u32)> = Vec::new();
    for i in start_codon_positions.iter() {
        let mut peptide = String::new();
        let orf_start: u32 = *i;
        let mut orf_end: u32 = *i;
        for codon_start in (*i..rna_sequence_.len() as u32).step_by(3) {
            if codon_start + 2 >= rna_sequence_.len() as u32 {
                break;
            }
            let codon: &str = &rna_sequence_[codon_start as usize..=(codon_start + 2) as usize].to_uppercase();
            let amino_acid: &str = CODON_TABLE[codon];
            peptide.push_str(amino_acid);
            orf_end = codon_start + 2;
            if amino_acid == "*" {
                break;
            }
        }
        peptides.push((peptide.clone().into_boxed_str(), orf_start, orf_end, peptide.len() as u32));
    }

    peptides
}

