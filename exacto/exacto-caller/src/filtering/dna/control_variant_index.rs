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


use exacto_core::prelude::{ReadID, ReadSupport, ReferenceChromosomeID, ReferencePosition};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::prelude::*;
use crate::calling::dna::variant_record_clustering::rebuild_orphan_pool_as_insertions;


pub struct DNAControlVariantIndex {
    /// HashMap<(chromosome ID, position 1), the bases the SNVs spell>
    position_snv_map: HashMap<(ReferenceChromosomeID, ReferencePosition), HashSet<String>>,

    /// HashMap<(chromosome ID, position 1 / bin size, index type), BTreeMap<position 1, graph operations>>
    position_1_map: HashMap<(ReferenceChromosomeID, u32, VariantType), BTreeMap<ReferencePosition, HashSet<Arc<GraphOperation>>>>,

    /// HashMap<(chromosome ID, position 2 / bin size, index type), BTreeMap<position 2, graph operations>>
    position_2_map: HashMap<(ReferenceChromosomeID, u32, VariantType), BTreeMap<ReferencePosition, HashSet<Arc<GraphOperation>>>>,

    bin_size: u32
}

impl DNAControlVariantIndex {
    pub fn new(
        variant_records: &[Arc<VariantRecord>],
        bin_size: u32,
        min_read_support: ReadSupport
    ) -> Self {
        // Step 1. Drop the records that too few reads carry. A read counts for the allele it
        // spells at its positions, whichever strand it was read from.
        let mut variant_records: Vec<Arc<VariantRecord>> = if min_read_support <= 1 {
            variant_records.to_vec()
        } else {
            let allele = |variant_record: &VariantRecord| -> (VariantType, ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, String) {
                (
                    variant_record.get_variant_type().clone(),
                    variant_record.get_chromosome_1(),
                    variant_record.get_position_1(),
                    variant_record.get_chromosome_2(),
                    variant_record.get_position_2(),
                    variant_record.get_standardized_sequence()
                )
            };
            let mut read_support: HashMap<(VariantType, ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, String), HashSet<ReadID>> = HashMap::new();
            for variant_record in variant_records.iter() {
                read_support
                    .entry(allele(variant_record))
                    .or_insert_with(HashSet::new)
                    .insert(variant_record.get_read_id());
            }
            variant_records
                .iter()
                .filter(|variant_record| {
                    read_support
                        .get(&allele(variant_record))
                        .map_or(false, |read_ids| read_ids.len() >= min_read_support as usize)
                })
                .cloned()
                .collect()
        };

        // Step 2. Give the clip records their insertion form as well.
        let clip_indices: Vec<usize> = (0..variant_records.len())
            .filter(|&i| variant_records[i].get_variant_type() == &VariantType::Breakpoint && !variant_records[i].is_resolved())
            .collect();
        let clip_insertions: Vec<Arc<VariantRecord>> = rebuild_orphan_pool_as_insertions(&variant_records, &clip_indices);
        variant_records.extend(clip_insertions);

        // Step 3. Index the SNVs by position, and the other records by position 1 and by position 2.
        let mut position_snv_map: HashMap<(ReferenceChromosomeID, ReferencePosition), HashSet<String>> = HashMap::new();
        let mut position_1_map: HashMap<(ReferenceChromosomeID, u32, VariantType), BTreeMap<ReferencePosition, HashSet<Arc<GraphOperation>>>> = HashMap::new();
        let mut position_2_map: HashMap<(ReferenceChromosomeID, u32, VariantType), BTreeMap<ReferencePosition, HashSet<Arc<GraphOperation>>>> = HashMap::new();
        for variant_record in variant_records.iter() {
            if variant_record.get_variant_type() == &VariantType::SingleNucleotideVariant {
                position_snv_map
                    .entry((variant_record.get_chromosome_1(), variant_record.get_position_1()))
                    .or_insert_with(HashSet::new)
                    .insert(variant_record.get_standardized_sequence());
                continue;
            }
            let graph_operation: Arc<GraphOperation> = Arc::new(variant_record.get_graph_operation().clone());
            let variant_type: VariantType = index_type(variant_record.get_variant_type());
            for (position_map, chromosome, position) in [
                (&mut position_1_map, variant_record.get_chromosome_1(), variant_record.get_position_1()),
                (&mut position_2_map, variant_record.get_chromosome_2(), variant_record.get_position_2())
            ] {
                position_map
                    .entry((chromosome, position / bin_size, variant_type.clone()))
                    .or_insert_with(BTreeMap::new)
                    .entry(position)
                    .or_insert_with(HashSet::new)
                    .insert(graph_operation.clone());
            }
        }

        Self {
            position_snv_map: position_snv_map,
            position_1_map: position_1_map,
            position_2_map: position_2_map,
            bin_size: bin_size
        }
    }

    pub fn get_snv_bases(&self, chromosome: ReferenceChromosomeID, position: ReferencePosition) -> Option<&HashSet<String>> {
        self.position_snv_map.get(&(chromosome, position))
    }
    
    pub fn get_nearby_graph_operations<'s>(
        &'s self,
        graph_operation: &GraphOperation,
        max_distance: u32
    ) -> impl Iterator<Item = &'s GraphOperation> + 's {
        let variant_type: VariantType = index_type(graph_operation.get_variant_type());
        let bin_size: u32 = self.bin_size;
        [
            (&self.position_1_map, graph_operation.get_chromosome_1(), graph_operation.get_position_1()),
            (&self.position_2_map, graph_operation.get_chromosome_2(), graph_operation.get_position_2())
        ]
            .into_iter()
            .flat_map(move |(position_map, chromosome, position)| {
                let bin: u32 = position / bin_size;
                let variant_type: VariantType = variant_type.clone();
                (bin.saturating_sub(1)..=(bin + 1))
                    .filter_map(move |bin| position_map.get(&(chromosome, bin, variant_type.clone())))
                    .flat_map(move |btree| btree.range(position.saturating_sub(max_distance)..=position + max_distance))
                    .flat_map(|(_, graph_operations)| graph_operations.iter().map(Arc::as_ref))
            })
    }
}


fn index_type(variant_type: &VariantType) -> VariantType {
    match variant_type {
        VariantType::Translocation => VariantType::Breakpoint,
        _ => variant_type.clone()
    }
}
