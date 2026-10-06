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


use exacto_caller::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;


pub(crate) struct ClusterDepths {
    /// HashMap<chromosome ID, (first bases, last bases)>
    exons: HashMap<u16, (Vec<u32>, Vec<u32>)>
}

impl ClusterDepths {
    pub(crate) fn new(transcript_models: &[Arc<TranscriptModel>]) -> Self {
        let mut exons: HashMap<u16, (Vec<u32>, Vec<u32>)> = HashMap::new();
        for transcript_model in transcript_models.iter() {
            let mut blocks: Vec<(u16, u32, u32)> = transcript_model
                .get_exons()
                .iter()
                .map(|exon| (exon.reference_chromosome_id, exon.reference_start, exon.reference_end))
                .collect();
            blocks.sort_unstable();
            let mut merged: Vec<(u16, u32, u32)> = Vec::with_capacity(blocks.len());
            for (chromosome_id, start, end) in blocks {
                match merged.last_mut() {
                    Some(last) if last.0 == chromosome_id && start <= last.2.saturating_add(1) => last.2 = last.2.max(end),
                    _ => merged.push((chromosome_id, start, end))
                }
            }
            for (chromosome_id, start, end) in merged {
                let (starts, ends): &mut (Vec<u32>, Vec<u32>) = exons.entry(chromosome_id).or_default();
                starts.push(start);
                ends.push(end);
            }
        }
        for (starts, ends) in exons.values_mut() {
            starts.sort_unstable();
            ends.sort_unstable();
        }
        Self { exons }
    }
    
    fn at(&self, chromosome_id: u16, position: u32) -> u32 {
        self.exons.get(&chromosome_id).map_or(0, |(starts, ends)| {
            (starts.partition_point(|&start| start <= position) - ends.partition_point(|&end| end < position)) as u32
        })
    }

    pub(crate) fn at_site(&self, op: &GraphOperation, carriers: &[&TranscriptModel]) -> (u32, u32) {
        let (chromosome_1, position_1): (u16, u32) = (op.get_chromosome_1(), op.get_position_1());
        let (chromosome_2, position_2): (u16, u32) = (op.get_chromosome_2(), op.get_position_2());
        let at = |chromosome_id: u16, position: u32| -> u32 {
            let num_carriers_elsewhere: usize = carriers
                .iter()
                .filter(|transcript_model| {
                    !transcript_model.get_exons().iter().any(|exon| {
                        exon.reference_chromosome_id == chromosome_id
                            && exon.reference_start <= position
                            && position <= exon.reference_end
                    })
                })
                .count();
            self.at(chromosome_id, position) + num_carriers_elsewhere as u32
        };
        match op.get_variant_type() {
            VariantType::SingleNucleotideVariant | VariantType::MultiNucleotideVariant => {
                let depth: u32 = (position_1.min(position_2) + 1..position_1.max(position_2))
                    .map(|position| at(chromosome_1, position))
                    .max()
                    .unwrap_or(0);
                (depth, depth)
            },
            VariantType::Insertion | VariantType::Deletion => {
                let depth: u32 = at(chromosome_1, position_1).max(at(chromosome_2, position_2));
                (depth, depth)
            },
            _ => (at(chromosome_1, position_1), at(chromosome_2, position_2))
        }
    }
}
