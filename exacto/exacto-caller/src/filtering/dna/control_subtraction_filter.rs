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


use edit_distance::edit_distance;
use exacto_core::prelude::{ReferenceChromosomeID, ReferencePosition};
use std::cmp::{min, max};

use crate::prelude::*;
use crate::calling::dna::variant_record_clustering::calculate_max_dna_clustering_distance;
use crate::filtering::dna::control_variant_index::DNAControlVariantIndex;


pub struct DNAControlSubtractionFilter<'a> {
    control_variant_index: &'a DNAControlVariantIndex,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64,
    apply_infinite_sites_assumption: bool,
    stranded: bool
}

impl <'a> DNAControlSubtractionFilter<'a> {
    pub fn new(
        control_variant_index: &'a DNAControlVariantIndex,
        min_size_proportion: f64,
        max_ins_norm_edit_distance: f64,
        max_clustering_distance: u32,
        sequencing_error: f64,
        apply_infinite_sites_assumption: bool,
        stranded: bool
    ) -> Self {
        Self {
            control_variant_index: control_variant_index,
            min_size_proportion: min_size_proportion,
            max_ins_norm_edit_distance: max_ins_norm_edit_distance,
            max_clustering_distance: max_clustering_distance,
            sequencing_error: sequencing_error,
            apply_infinite_sites_assumption: apply_infinite_sites_assumption,
            stranded: stranded
        }
    }
}

impl VariantFilter for DNAControlSubtractionFilter<'_> {
    type Input = VariantRecord;

    fn passes(&self, variant_record: &VariantRecord) -> bool {
        let graph_operation: &GraphOperation = variant_record.get_graph_operation();
        if graph_operation.get_variant_type() == &VariantType::SingleNucleotideVariant {
            // Under the infinite sites assumption a position mutates once, so any control SNV
            // there is this one; without it only the same base is.
            return !self.control_variant_index
                .get_snv_bases(graph_operation.get_chromosome_1(), graph_operation.get_position_1())
                .is_some_and(|bases| {
                    self.apply_infinite_sites_assumption || bases.contains(&graph_operation.get_standardized_sequence())
                });
        }
        self.control_variant_index
            .get_nearby_graph_operations(graph_operation, self.max_clustering_distance)
            .all(|control_graph_operation| self.represent_different_dna_variants(graph_operation, control_graph_operation))
    }
}

impl DNAControlSubtractionFilter<'_> {
    fn represent_different_dna_variants(&self, a: &GraphOperation, b: &GraphOperation) -> bool {
        // Fast early exit for different chromosomes.
        if (a.get_chromosome_1() != b.get_chromosome_1() && a.get_chromosome_2() != b.get_chromosome_2()) &&
            (a.get_chromosome_1() != b.get_chromosome_2() && a.get_chromosome_2() != b.get_chromosome_1()) {
            return true;
        }

        // Strand comparison.
        if self.stranded {
            if a.get_strand_1() != b.get_strand_1() || a.get_strand_2() != b.get_strand_2() {
                return true;
            }
        }

        // A breakpoint and a translocation share an index, and are two variants unless one of them
        // is a clip (operation 2 Noop): a clip is a Breakpoint record that may be one side of a
        // translocation.
        let is_clip = |operation: &GraphOperation| *operation.get_operation_type_2() == GraphOperationType::Noop;
        if a.get_variant_type() != b.get_variant_type() && is_clip(a) == is_clip(b) {
            return true;
        }

        // Precompute reusable properties.
        let size_a: isize = a.get_variant_size();
        let size_b: isize = b.get_variant_size();
        let max_size: isize = max(size_a, size_b);
        let min_size: isize = min(size_a, size_b);
        let size_proportion: f64 = if max_size > 0 {
            min_size as f64 / max_size as f64
        } else {
            1.0
        };

        // Infinite sites assumption.
        if self.apply_infinite_sites_assumption {
            if (a.get_chromosome_1() == b.get_chromosome_1() && a.get_position_1() == b.get_position_1()) ||
                (a.get_chromosome_2() == b.get_chromosome_2() && a.get_position_2() == b.get_position_2()) {
                return false;
            }
        }

        // A clip is a breakpoint seen from one side, and it is the junction it clips at when its one
        // breakend is one of the junction's, by chromosome and operation.
        if is_clip(a) != is_clip(b) {
            let (clip, junction): (&GraphOperation, &GraphOperation) = if is_clip(a) { (a, b) } else { (b, a) };
            let max_distance: u32 = calculate_max_dna_clustering_distance(
                clip.get_sequence_length() as u32,
                self.sequencing_error,
                self.max_clustering_distance
            );
            let meets = |chromosome: ReferenceChromosomeID, position: ReferencePosition, operation: &GraphOperationType| -> bool {
                clip.get_chromosome_1() == chromosome
                    && clip.get_operation_type_1() == operation
                    && clip.get_position_1().abs_diff(position) <= max_distance
            };
            return !meets(junction.get_chromosome_1(), junction.get_position_1(), junction.get_operation_type_1())
                && !meets(junction.get_chromosome_2(), junction.get_position_2(), junction.get_operation_type_2());
        }

        let pos1_distance: u32 = a.get_position_1().abs_diff(b.get_position_1());
        let pos2_distance: u32 = a.get_position_2().abs_diff(b.get_position_2());

        // Translocation.
        if a.get_variant_type().clone() == VariantType::Translocation {
            let cross_pos1_distance: u32 = a.get_position_1().abs_diff(b.get_position_2());
            let cross_pos2_distance: u32 = a.get_position_2().abs_diff(b.get_position_1());
            if (a.get_chromosome_1() == b.get_chromosome_1() && pos1_distance <= self.max_clustering_distance) ||
                (a.get_chromosome_2() == b.get_chromosome_2() && pos2_distance <= self.max_clustering_distance) ||
                (a.get_chromosome_1() == b.get_chromosome_2() && cross_pos1_distance <= self.max_clustering_distance) ||
                (a.get_chromosome_2() == b.get_chromosome_1() && cross_pos2_distance <= self.max_clustering_distance) {
                return false;
            } else {
                return true;
            }
        }

        let max_distance: u32 = calculate_max_dna_clustering_distance(
            max(size_a, size_b) as u32,
            self.sequencing_error,
            self.max_clustering_distance
        );

        // Insertion.
        if a.get_variant_type().clone() == VariantType::Insertion {
            let a_sequence: String = a.get_standardized_sequence();
            let b_sequence: String = b.get_standardized_sequence();
            let edit_distance: f64 = edit_distance(a_sequence.as_str(), b_sequence.as_str()) as f64;
            let max_seq_length: f64 = max(a.get_sequence_length(), b.get_sequence_length()) as f64;
            let normalized_edit_distance: f64 = edit_distance / max_seq_length;
            if a.get_chromosome_1() == b.get_chromosome_1() &&
                pos1_distance <= max_distance &&
                size_proportion >= self.min_size_proportion &&
                normalized_edit_distance <= self.max_ins_norm_edit_distance {
                return false;
            } else {
                return true;
            }
        }

        // An MNV that does not share a position under the infinite sites assumption is the same
        // variant only when it spells the same bases.
        if a.get_variant_type() == &VariantType::MultiNucleotideVariant
            && a.get_standardized_sequence() != b.get_standardized_sequence() {
            return true;
        }

        // Deletion, breakpoint or MNV.
        if (a.get_chromosome_1() == b.get_chromosome_1() && pos1_distance <= max_distance && size_proportion >= self.min_size_proportion) ||
            (a.get_chromosome_2() == b.get_chromosome_2() && pos2_distance <= max_distance && size_proportion >= self.min_size_proportion) {
            return false;
        }

        true
    }
}


#[cfg(test)]
#[path = "../../tests/filtering/dna/control_subtraction_filter.rs"]
mod tests;
