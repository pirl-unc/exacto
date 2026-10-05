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
use exacto_caller::prelude::*;
use exacto_core::log_info;
use std::collections::HashSet;

use crate::io::records::RNAReadClusterVariantRecord;


/// A DNA call on the BAM's chromosome ids:
/// (chromosome 1, chromosome 2, position 1, position 2, operation 1, operation 2, sequence).
type DNAVariantKey<'a> = (u16, u16, u32, u32, GraphOperationType, GraphOperationType, &'a str);

/// An allowed variant on the BAM's chromosome ids:
/// (chromosome 1, chromosome 2, position 1, position 2, variant type, sequence).
type AllowedVariant = (u16, u16, u32, u32, VariantType, Box<str>);


pub(crate) struct KnownVariants<'a> {
    dna_variants: HashSet<DNAVariantKey<'a>>,
    allowed_variants: Option<Vec<AllowedVariant>>
}

impl<'a> KnownVariants<'a> {
    pub(crate) fn new(
        dna_variant_records: Option<&'a Vec<DNAVariantRecord>>,
        allowed_variant_records: Option<&Vec<RNAReadClusterVariantRecord>>,
        chromosome_names_map: &BiMap<Box<str>, u16>
    ) -> Self {
        let mut dna_variants: HashSet<DNAVariantKey<'a>> = HashSet::new();
        if let Some(records) = dna_variant_records {
            log_info!("Creating a HashSet of the DNA variant records.");
            dna_variants.reserve(records.len());
            let mut num_absent: usize = 0;
            for record in records {
                let (Some(&chromosome_1), Some(&chromosome_2)) = (
                    chromosome_names_map.get_by_left(record.chromosome_1.as_ref()),
                    chromosome_names_map.get_by_left(record.chromosome_2.as_ref())
                ) else {
                    num_absent += 1;
                    continue;
                };
                dna_variants.insert((
                    chromosome_1,
                    chromosome_2,
                    record.position_1,
                    record.position_2,
                    record.operation_1
                        .parse::<GraphOperationType>()
                        .expect("Invalid DNA operation_1"),
                    record.operation_2
                        .parse::<GraphOperationType>()
                        .expect("Invalid DNA operation_2"),
                    record.sequence.as_ref(),
                ));
            }
            if num_absent > 0 {
                log_info!(
                    "Left out {} of {} DNA variant record(s) on a contig the RNA BAM file does not hold.",
                    num_absent,
                    records.len()
                );
            }
        }

        let allowed_variants: Option<Vec<AllowedVariant>> = allowed_variant_records.map(|records| {
            let mut num_absent: usize = 0;
            let allowed_variants: Vec<AllowedVariant> = records
                .iter()
                .filter_map(|record| {
                    let (Some(&chromosome_1), Some(&chromosome_2)) = (
                        chromosome_names_map.get_by_left(record.chromosome_1.as_ref()),
                        chromosome_names_map.get_by_left(record.chromosome_2.as_ref())
                    ) else {
                        num_absent += 1;
                        return None;
                    };
                    let variant_type: VariantType = record.variant_type
                        .parse::<VariantType>()
                        .unwrap_or_else(|_| panic!("Invalid variant_type of an allowed variant: {}", record.variant_type));
                    Some((chromosome_1, chromosome_2, record.position_1, record.position_2, variant_type, record.sequence.clone()))
                })
                .collect();
            if num_absent > 0 {
                log_info!(
                    "Left out {} of {} allowed variant(s) on a contig the RNA BAM file does not hold.",
                    num_absent,
                    records.len()
                );
            }
            log_info!("Restricting the variant calls to {} allowed variant(s).", allowed_variants.len());
            allowed_variants
        });

        Self { dna_variants, allowed_variants }
    }
    
    pub(crate) fn confirms(&self, op: &GraphOperation) -> bool {
        !self.dna_variants.is_empty() && self.dna_variants.contains(&(
            op.get_chromosome_1(),
            op.get_chromosome_2(),
            op.get_position_1(),
            op.get_position_2(),
            op.get_operation_type_1().clone(),
            op.get_operation_type_2().clone(),
            op.get_standardized_sequence().as_str()
        ))
    }
    
    pub(crate) fn allows(&self, op: &GraphOperation, max_distance: u32, max_ins_norm_edit_distance: f64) -> bool {
        let Some(allowed_variants) = self.allowed_variants.as_ref() else {
            return true;
        };
        allowed_variants.iter().any(|(chromosome_1, chromosome_2, position_1, position_2, variant_type, sequence)| {
            variant_type == op.get_variant_type()
                && *chromosome_1 == op.get_chromosome_1()
                && *chromosome_2 == op.get_chromosome_2()
                && position_1.abs_diff(op.get_position_1()) <= max_distance
                && position_2.abs_diff(op.get_position_2()) <= max_distance
                && (*variant_type != VariantType::Insertion || {
                    let num_edits: usize = edit_distance::edit_distance(sequence, op.get_sequence());
                    num_edits as f64 <= max_ins_norm_edit_distance * sequence.len().max(op.get_sequence().len()) as f64
                })
        })
    }
}


#[cfg(test)]
#[path = "../tests/variant_filtering/known_variants.rs"]
mod tests;