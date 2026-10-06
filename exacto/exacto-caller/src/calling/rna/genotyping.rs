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


use exacto_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;


use crate::prelude::*;


pub fn genotype_rna_variants(
    read_ids: &HashSet<ReadID>,
    variant_calls: &[VariantCall],
    transcript_models_map: &HashMap<ReadID, Arc<TranscriptModel>>
) -> HashMap<VariantID, Vec<(ReadID, Allele)>> {
    // Step 1. Per read, the two things the ladder asks of it: its aligned intervals, and
    // the left flank of every insertion it carries. A read without a model has neither
    // and lands on NotCovered at every site.
    let mut read_exons: HashMap<ReadID, Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)>> = HashMap::new();
    let mut read_insertion_anchors: HashMap<ReadID, HashSet<(ReferenceChromosomeID, ReferencePosition)>> = HashMap::new();
    for &read_id in read_ids.iter() {
        let Some(transcript_model) = transcript_models_map.get(&read_id) else {
            continue;
        };
        read_exons.insert(
            read_id,
            transcript_model
                .get_exons()
                .iter()
                .map(|exon| (exon.reference_chromosome_id, exon.reference_start, exon.reference_end))
                .collect()
        );
        let anchors: HashSet<(ReferenceChromosomeID, ReferencePosition)> = transcript_model
            .get_variant_records()
            .iter()
            .filter(|record| *record.get_variant_type() == VariantType::Insertion)
            .map(|record| (record.get_chromosome_1(), record.get_position_1().min(record.get_position_2())))
            .collect();
        if !anchors.is_empty() {
            read_insertion_anchors.insert(read_id, anchors);
        }
    }

    // Step 2. One allele per read per call.
    // HashMap<variant call ID, Vec<(read ID, Allele)>>
    let mut genotypes: HashMap<VariantID, Vec<(ReadID, Allele)>> = HashMap::with_capacity(variant_calls.len());
    for variant_call in variant_calls.iter() {
        let site: SiteSpan = SiteSpan::from_variant_call(variant_call);
        let anchor: (ReferenceChromosomeID, ReferencePosition) = (site.chromosome_id, site.reference_start);
        let alternate_read_ids: HashSet<ReadID> = variant_call.get_read_ids().into_iter().collect();

        let alleles: Vec<(ReadID, Allele)> = read_ids
            .iter()
            .map(|&read_id| {
                let allele: Allele = if alternate_read_ids.contains(&read_id) {
                    Allele::Alternate
                } else if site.insertion_length.is_some()
                    && read_insertion_anchors.get(&read_id).is_some_and(|anchors| anchors.contains(&anchor))
                {
                    Allele::NotCovered
                } else if read_exons.get(&read_id).is_some_and(|exons| site.is_covered_by(exons)) {
                    Allele::Reference
                } else {
                    Allele::NotCovered
                };
                (read_id, allele)
            })
            .collect();
        genotypes.insert(variant_call.get_id(), alleles);
    }

    genotypes
}


struct SiteSpan {
    chromosome_id: ReferenceChromosomeID,
    reference_start: ReferencePosition,
    reference_end: ReferencePosition,
    insertion_length: Option<u32>,
    side_2: Option<(ReferenceChromosomeID, ReferencePosition)>
}

impl SiteSpan {
    fn from_variant_call(variant_call: &VariantCall) -> Self {
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        let (position_1, position_2): (ReferencePosition, ReferencePosition) = (operation.get_position_1(), operation.get_position_2());
        let (mut reference_start, mut reference_end): (ReferencePosition, ReferencePosition) = (position_1.min(position_2), position_1.max(position_2));
        if matches!(
            operation.get_variant_type(),
            VariantType::SingleNucleotideVariant | VariantType::MultiNucleotideVariant
        ) && reference_end - reference_start >= 2 {
            reference_start += 1;
            reference_end -= 1;
        }
        let is_breakend: bool = matches!(
            operation.get_variant_type(),
            VariantType::Breakpoint | VariantType::Translocation | VariantType::FusionGene | VariantType::CircularRNA
        );
        if is_breakend {
            (reference_start, reference_end) = (position_1, position_1);
        }
        Self {
            chromosome_id: operation.get_chromosome_1(),
            reference_start,
            reference_end,
            insertion_length: (*operation.get_variant_type() == VariantType::Insertion)
                .then(|| operation.get_sequence_length() as u32),
            side_2: is_breakend.then(|| (operation.get_chromosome_2(), position_2))
        }
    }

    fn is_covered_by(&self, exons: &[(ReferenceChromosomeID, ReferencePosition, ReferencePosition)]) -> bool {
        let covers_side_1: bool = exons.iter().any(|&(chromosome_id, start, end)| {
            chromosome_id == self.chromosome_id
                && start <= self.reference_start
                && match self.insertion_length {
                Some(insertion_length) => self.reference_start <= end.saturating_add(insertion_length),
                None => self.reference_end <= end
            }
        });
        covers_side_1 || self.side_2.is_some_and(|(chromosome_id_2, position_2)| {
            exons.iter().any(|&(chromosome_id, start, end)| {
                chromosome_id == chromosome_id_2 && start <= position_2 && position_2 <= end
            })
        })
    }
}


#[cfg(test)]
#[path = "../../tests/calling/rna/genotyping.rs"]
mod tests;