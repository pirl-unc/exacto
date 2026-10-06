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
use exacto_core::prelude::{FastaMap, Strand};
use std::collections::HashMap;


pub(crate) struct PlacedJunction {
    pub(crate) read_position: u32,
    pub(crate) number: u16,
    pub(crate) junction: SpliceJunction,
    windows: [Option<(u32, u32)>; 2]
}


pub(crate) fn place_read_junctions(
    transcript_model: &TranscriptModel,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> Vec<PlacedJunction> {
    // (low, high, sequence) of every insertion the read carries.
    let insertions: Vec<(u32, u32, &str)> = transcript_model
        .get_variant_records()
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::Insertion)
        .map(|variant_record| (
            variant_record.get_position_1().min(variant_record.get_position_2()),
            variant_record.get_position_1().max(variant_record.get_position_2()),
            variant_record.get_sequence()
        ))
        .collect();

    transcript_model
        .get_splice_junctions()
        .iter()
        .map(|junction| {
            let position_1: u32 = junction.reference_position_1;
            let position_2: u32 = junction.reference_position_2;
            let forward: bool = position_1 < position_2;
            let chromosome: Option<&Box<str>> = chromosome_names_map.get_by_right(&junction.reference_chromosome_id_1);

            // The insertion that frees a boundary sits between that boundary and the exon base
            // next to it, so it is found by the pair of reference positions it lies between.
            let mut placed: [u32; 2] = [position_1, position_2];
            let mut windows: [Option<(u32, u32)>; 2] = [None, None];
            for (index, boundary) in [JunctionBoundary::Donor, JunctionBoundary::Acceptor].into_iter().enumerate() {
                let exon_is_below: bool = match boundary {
                    JunctionBoundary::Donor => forward,
                    JunctionBoundary::Acceptor => !forward
                };
                let flanks: (u32, u32) = if exon_is_below {
                    (placed[index].saturating_sub(1), placed[index])
                } else {
                    (placed[index], placed[index].saturating_add(1))
                };
                let (Some(chromosome), Some((_, _, sequence))) = (
                    chromosome,
                    insertions.iter().find(|(low, high, _)| (*low, *high) == flanks)
                ) else {
                    continue;
                };
                let slid: u32 = normalise_junction_boundary(position_1, position_2, sequence, boundary, chromosome, fasta_map);
                windows[index] = Some((placed[index].min(slid), placed[index].max(slid)));
                placed[index] = slid;
            }

            PlacedJunction {
                read_position: junction.read_position_1,
                number: junction.splice_junction_number,
                junction: SpliceJunction::new(
                    junction.reference_chromosome_id_1,
                    junction.reference_chromosome_id_2,
                    placed[0],
                    placed[1],
                    junction.reference_strand_1.clone(),
                    junction.reference_strand_2.clone()
                ),
                windows
            }
        })
        .collect()
}


pub(crate) fn pool_insertion_abutting_junctions(reads: &mut [Vec<PlacedJunction>]) -> usize {
    // A boundary's group is everything about its junction except the boundary itself.
    type GroupKey = (u16, u16, u32, Strand, Strand);

    let mut num_moved: usize = 0;
    for (index, boundary) in [JunctionBoundary::Donor, JunctionBoundary::Acceptor].into_iter().enumerate() {
        let group_key = |junction: &SpliceJunction| -> GroupKey {
            let opposite: u32 = match boundary {
                JunctionBoundary::Donor => junction.position_2,
                JunctionBoundary::Acceptor => junction.position_1
            };
            (junction.chromosome_1, junction.chromosome_2, opposite, junction.strand_1.clone(), junction.strand_2.clone())
        };
        let placement = |junction: &SpliceJunction| -> u32 {
            match boundary {
                JunctionBoundary::Donor => junction.position_1,
                JunctionBoundary::Acceptor => junction.position_2
            }
        };

        // Step 1. Tally: one vote per junction at its placement, and the union of the slide
        // windows seen there.
        // HashMap<group, HashMap<placement, ((window low, window high), votes)>>
        let mut groups: HashMap<GroupKey, HashMap<u32, ((u32, u32), usize)>> = HashMap::new();
        for placed in reads.iter().flatten() {
            let at: u32 = placement(&placed.junction);
            let window: (u32, u32) = placed.windows[index].unwrap_or((at, at));
            let entry: &mut ((u32, u32), usize) = groups
                .entry(group_key(&placed.junction))
                .or_default()
                .entry(at)
                .or_insert((window, 0));
            entry.0.0 = entry.0.0.min(window.0);
            entry.0.1 = entry.0.1.max(window.1);
            entry.1 += 1;
        }

        // Step 2. Per group, single-link the placements whose windows reach one another, then
        // send every component to its mode. Ties go to the larger placement, so the winner
        // never depends on iteration order.
        // HashMap<(group, losing placement), modal placement>
        let mut remap: HashMap<(GroupKey, u32), u32> = HashMap::new();
        for (group, placements) in groups.into_iter() {
            let mut ordered: Vec<(u32, (u32, u32), usize)> = placements
                .into_iter()
                .map(|(at, (window, votes))| (at, window, votes))
                .collect();
            ordered.sort_unstable_by_key(|entry| entry.0);
            let mut components: Vec<Vec<(u32, (u32, u32), usize)>> = Vec::new();
            let mut reach: u32 = 0;
            for entry in ordered.into_iter() {
                let (_, (window_low, window_high), _) = entry;
                if components.is_empty() || window_low > reach {
                    components.push(vec![entry]);
                    reach = window_high;
                } else {
                    components.last_mut().unwrap().push(entry);
                    reach = reach.max(window_high);
                }
            }
            for component in components.into_iter() {
                let modal: u32 = component
                    .iter()
                    .max_by_key(|(at, _, votes)| (*votes, *at))
                    .map(|(at, _, _)| *at)
                    .unwrap();
                for (at, _, _) in component.into_iter() {
                    if at != modal {
                        remap.insert((group.clone(), at), modal);
                    }
                }
            }
        }

        // Step 3. Move. Only a boundary that had an insertion beside it can move; one with no
        // window voted in the tally but stays where it is.
        for placed in reads.iter_mut().flatten() {
            if placed.windows[index].is_none() {
                continue;
            }
            let key: (GroupKey, u32) = (group_key(&placed.junction), placement(&placed.junction));
            if let Some(&modal) = remap.get(&key) {
                match boundary {
                    JunctionBoundary::Donor => placed.junction.position_1 = modal,
                    JunctionBoundary::Acceptor => placed.junction.position_2 = modal
                }
                num_moved += 1;
            }
        }
    }
    num_moved
}
