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


use bimap::BiMap;
use exacto_core::prelude::*;

use crate::prelude::{GraphOperation, VariantType};


pub fn is_repeat_variant(
    graph_operation: &GraphOperation,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    min_homopolymer_len: u32,
    min_dinucleotide_context_len: u32
) -> (bool, u32) {
    let variant_type: &VariantType = graph_operation.get_variant_type();
    if !matches!(
        variant_type,
        VariantType::Insertion |
        VariantType::Deletion |
        VariantType::SingleNucleotideVariant |
        VariantType::MultiNucleotideVariant
    ) {
        return (false, 0);
    }

    let chromosome: &str = chromosome_names_map
        .get_by_right(&graph_operation.get_chromosome_1())
        .unwrap();

    let position_1: ReferencePosition = graph_operation.get_position_1();
    let position_2: ReferencePosition = graph_operation.get_position_2();

    let get_fasta_sequence = |start: ReferencePosition, end: ReferencePosition| -> Box<str> {
        if start >= 1 && end <= fasta_map.get_length(chromosome) as u32 {
            fasta_map
                .get_sequence(chromosome, start as usize, end as usize)
                .to_uppercase()
                .into_boxed_str()
        } else {
            "".into()
        }
    };

    // Count the number of positions that matches the given single-nucleotide
    // motif to the left of position (towards 5') in the reference genome.
    let count_left_homopolymer = |mut position: ReferencePosition, base: &str| -> u32 {
        let mut n: u32 = 0;
        while position >= 1 {
            if &*get_fasta_sequence(position, position) == base {
                n += 1;
                if position == 1 {
                    break;
                }
                position -= 1;
            } else {
                break;
            }
        }
        n
    };

    // Count the number of positions that matches the given single-nucleotide
    // motif to the right of position (towards 3') in the reference genome.
    let count_right_homopolymer = |mut position: ReferencePosition, base: &str| -> u32 {
        let mut n: u32 = 0;
        loop {
            if &*get_fasta_sequence(position, position) == base {
                n += 1;
                position = position.saturating_add(1);
            } else {
                break;
            }
        }
        n
    };

    // Count the number of positions that matches the given dinucleotide
    // motif to the left of position (towards 5') in the reference genome.
    let count_left_dinucleotide_repeat = |mut position: ReferencePosition, motif: &str| -> u32 {
        let mut n: u32 = 0;
        while position >= 2 {
            let s: Box<str> = get_fasta_sequence(position - 1, position);
            if s.len() != 2 {
                break;
            }
            if &*s == motif {
                n += 2;
                position -= 2;
            } else {
                break;
            }
        }
        n
    };

    // Count the number of positions that matches the given dinucleotide
    // motif to the right of position (towards 3') in the reference genome.
    let count_right_dinucleotide_repeat = |mut position: ReferencePosition, motif: &str| -> u32 {
        let mut n: u32 = 0;
        loop {
            let s: Box<str> = get_fasta_sequence(position, position + 1);
            if s.len() != 2 {
                break;
            }
            if &*s == motif {
                n += 2;
                position = position.saturating_add(2);
            } else {
                break;
            }
        }
        n
    };

    // Length of the longest homopolymer run in the reference genome that overlaps the
    // inclusive window `[start, end]`.
    let homopolymer_context = |start: ReferencePosition, end: ReferencePosition| -> u32 {
        let mut longest: u32 = 0;
        for position in start..=end {
            let base: Box<str> = get_fasta_sequence(position, position);
            if base.is_empty() {
                continue;
            }
            // The run through `position`: the base itself plus its matching neighbours.
            let length: u32 = 1
                + count_left_homopolymer(position.saturating_sub(1), &*base)
                + count_right_homopolymer(position + 1, &*base);
            longest = longest.max(length);
        }
        longest
    };

    // Insertion.
    if variant_type == &VariantType::Insertion {
        let sequence: String = graph_operation.get_standardized_sequence();
        assert!(sequence.len() > 0);

        // Homopolymer insertion.
        if is_homopolymer_sequence(sequence.as_str()) {
            let base: &str = sequence.chars().next().map(|c| &sequence[..c.len_utf8()]).unwrap();

            // Count the length of the left homopolymer sequence, if any.
            let left_homopolymer_length: u32 = count_left_homopolymer(position_1, base);

            // Count the length of the right homopolymer sequence, if any.
            let right_homopolymer_length: u32 = count_right_homopolymer(position_2, base);

            // The run the insertion leaves is held to the minimum the deletion branch holds the
            // run a deletion takes bases out of: inserting `G` beside a single `G` makes a run
            // of 2, which is not a repeat.
            let homopolymer_length: u32 = left_homopolymer_length + sequence.len() as u32 + right_homopolymer_length;
            if (left_homopolymer_length > 0 || right_homopolymer_length > 0) && homopolymer_length >= min_homopolymer_len {
                return (true, homopolymer_length);
            }
        }

        // Dinucleotide repeat insertion.
        if sequence.len() == 1 {
            // Count the length of the left dinucleotide repeat sequence, if any.
            let left_sequence: Box<str> = get_fasta_sequence(position_1.saturating_sub(1), position_1);
            let left_repeat_length: u32 = count_left_dinucleotide_repeat(position_1, &*left_sequence);

            // Count the length of the right dinucleotide repeat sequence, if any.
            let right_sequence: Box<str> = get_fasta_sequence(position_2, position_2 + 1);
            let right_repeat_length: u32 = count_right_dinucleotide_repeat(position_2, &*right_sequence);

            // Both counters are anchored on the flanking motif itself, so each returns at
            // least 2 unconditionally; only a run that carries past the anchor is a repeat.
            if left_repeat_length >= min_dinucleotide_context_len ||
                right_repeat_length >= min_dinucleotide_context_len {
                return (true, left_repeat_length + sequence.len() as u32 + right_repeat_length);
            }
        } else {
            let result: Option<String> = get_dinucleotide_repeat_motif(sequence.as_str());
            if result.is_some() {
                let motif: String = result.unwrap();
                let motif_reversed: String = motif.chars().rev().collect();

                // Count the length of the left dinucleotide repeat sequence, if any.
                let left_repeat_length_1: u32 = count_left_dinucleotide_repeat(position_1, motif.as_str());
                let left_repeat_length_2: u32 = count_left_dinucleotide_repeat(position_1, motif_reversed.as_str());

                // Count the length of the right dinucleotide repeat sequence, if any.
                let right_repeat_length_1: u32 = count_right_dinucleotide_repeat(position_2, motif.as_str());
                let right_repeat_length_2: u32 = count_right_dinucleotide_repeat(position_2, motif_reversed.as_str());

                let length_1: u32 = if left_repeat_length_1 > 0 || right_repeat_length_1 > 0 {
                    left_repeat_length_1 + sequence.len() as u32 + right_repeat_length_1
                } else {
                    0
                };

                let length_2: u32 = if left_repeat_length_2 > 0 || right_repeat_length_2 > 0 {
                    left_repeat_length_2 + sequence.len() as u32 + right_repeat_length_2
                } else {
                    0
                };

                if length_1 > 0 || length_2 > 0 {
                    return (true, length_1.max(length_2));
                }
            }
        }
    }

    // Deletion.
    if variant_type == &VariantType::Deletion {
        let sequence: Box<str> = get_fasta_sequence(position_1 + 1, position_2 - 1);

        assert!(
            !sequence.is_empty(),
            "deletion at {}:{}-{} spans no reference bases",
            chromosome, position_1, position_2
        );

        // Homopolymer deletion.
        let deleted_length: u32 = sequence.len() as u32;
        let homopolymer_length: u32 = homopolymer_context(position_1 + 1, position_1 + 1)
            .max(homopolymer_context(position_2 - 1, position_2 - 1));
        if homopolymer_length >= min_homopolymer_len && deleted_length < homopolymer_length {
            return (true, homopolymer_length);
        }

        // Dinucleotide repeat deletion.
        if sequence.len() == 1 {
            // Count the length of the left dinucleotide repeat sequence, if any.
            let left_sequence: Box<str> = get_fasta_sequence(position_1.saturating_sub(1), position_1);
            let left_repeat_length: u32 = count_left_dinucleotide_repeat(position_1, &*left_sequence);

            // Count the length of the right dinucleotide repeat sequence, if any.
            let right_sequence: Box<str> = get_fasta_sequence(position_2, position_2 + 1);
            let right_repeat_length: u32 = count_right_dinucleotide_repeat(position_2, &*right_sequence);

            // Anchored on the flanking motif, as in the insertion branch: 2 is no evidence.
            if left_repeat_length >= min_dinucleotide_context_len ||
                right_repeat_length >= min_dinucleotide_context_len {
                return (true, left_repeat_length + sequence.len() as u32 + right_repeat_length);
            }
        } else {
            let result: Option<String> = get_dinucleotide_repeat_motif(&*sequence);
            if result.is_some() {
                let motif: String = result.unwrap();
                let motif_reversed: String = motif.chars().rev().collect();

                // Count the length of the left dinucleotide repeat sequence, if any.
                let left_repeat_length_1: u32 = count_left_dinucleotide_repeat(position_1, motif.as_str());
                let left_repeat_length_2: u32 = count_left_dinucleotide_repeat(position_1, motif_reversed.as_str());

                // Count the length of the right dinucleotide repeat sequence, if any.
                let right_repeat_length_1: u32 = count_right_dinucleotide_repeat(position_2, motif.as_str());
                let right_repeat_length_2: u32 = count_right_dinucleotide_repeat(position_2, motif_reversed.as_str());

                let length_1: u32 = if left_repeat_length_1 > 0 || right_repeat_length_1 > 0 {
                    left_repeat_length_1 + sequence.len() as u32 + right_repeat_length_1
                } else {
                    0
                };

                let length_2: u32 = if left_repeat_length_2 > 0 || right_repeat_length_2 > 0 {
                    left_repeat_length_2 + sequence.len() as u32 + right_repeat_length_2
                } else {
                    0
                };

                if length_1 > 0 || length_2 > 0 {
                    return (true, length_1.max(length_2));
                }
            }
        }
    }

    // Single and multi-nucleotide variants.
    if variant_type == &VariantType::SingleNucleotideVariant ||
        variant_type == &VariantType::MultiNucleotideVariant {
        let sequence: Box<str> = get_fasta_sequence(position_1 + 1, position_2 - 1);

        // Check if the substituted bases lie in a homopolymer repeat context.
        let homopolymer_length: u32 = homopolymer_context(position_1 + 1, position_2 - 1);
        if homopolymer_length >= min_homopolymer_len {
            return (true, homopolymer_length);
        }

        // Check if the substituted bases lie in a dinucleotide repeat context.
        let left_sequence: Box<str> = get_fasta_sequence(position_1.saturating_sub(1), position_1);
        let left_repeat_length: u32 = count_left_dinucleotide_repeat(position_1, &*left_sequence);

        let right_sequence: Box<str> = get_fasta_sequence(position_2, position_2 + 1);
        let right_repeat_length: u32 = count_right_dinucleotide_repeat(position_2, &*right_sequence);

        if left_repeat_length >= min_dinucleotide_context_len ||
            right_repeat_length >= min_dinucleotide_context_len {
            return (true, left_repeat_length + sequence.len() as u32 + right_repeat_length);
        }
    }

    (false, 0)
}


#[cfg(test)]
#[path = "../tests/filtering/repeat_context.rs"]
mod tests;
