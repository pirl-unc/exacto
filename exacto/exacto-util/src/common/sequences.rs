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
pub fn find_kmers(sequence: &str, k: usize) -> HashMap<Box<str>,Vec<usize>> {
    let mut kmer_index: HashMap<Box<str>,Vec<usize>> = HashMap::new();
    for i in 0..=sequence.len() - k {
        let kmer = &sequence[i..i + k];
        kmer_index.entry(kmer.to_string().into_boxed_str())
            .or_insert(Vec::new())
            .push(i);
    }
    kmer_index
}

/// Find substring start positions.
pub fn find_substring_positions(text: &str, substring: &str) -> Vec<usize> {
    let mut start_positions: Vec<usize> = Vec::new();
    let substring_len: usize = substring.len();
    if substring_len == 0 || substring_len > text.len() {
        return start_positions;
    }
    for (index, window) in text.as_bytes().windows(substring_len).enumerate() {
        if window == substring.as_bytes() {
            start_positions.push(index);
        }
    }
    start_positions
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
    assert!(is_valid_nucleotide_sequence(sequence), "The nucleotide sequence is not valid: {}.", sequence);
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

/// Translate a RNA sequence to all possible ORF peptides.
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
) -> Vec<(Box<str>,usize,usize,usize)> {
    let rna_sequence_: Box<str> = rna_sequence.replace('T', "U").replace('t', "u").into_boxed_str();
    let mut peptides: Vec<(Box<str>,usize,usize,usize)> = Vec::new();
    for frame in 0..3 {
        let mut peptide = String::new();
        let mut in_orf = false;
        let mut orf_start: usize = 0;
        let mut orf_end: usize = 0;
        for codon_start in (frame..rna_sequence.len()).step_by(3) {
            if codon_start + 3 > rna_sequence_.len() {
                break;
            }
            let codon: &str = &rna_sequence_[codon_start..codon_start + 3].to_uppercase();
            if start_codons.contains(codon) && !in_orf {
                // Start of a new ORF
                peptide = String::new();
                in_orf = true;
                orf_start = codon_start;
            }
            if in_orf {
                let amino_acid: &str = CODON_TABLE[codon];
                peptide.push_str(amino_acid);
                orf_end = codon_start + 2;
                if amino_acid == "*" {
                    if !peptide.is_empty() {
                        peptides.push((peptide.clone().into_boxed_str(), orf_start, orf_end, peptide.len()));
                    }
                    in_orf = false;
                }
            }
        }

        // If still in ORF after processing the sequence
        if in_orf && !peptide.is_empty() {
            peptides.push((peptide.clone().into_boxed_str(), orf_start, orf_end, peptide.len()));
        }
    }

    peptides
}

