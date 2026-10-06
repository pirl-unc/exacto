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
use exacto_core::prelude::{reverse_complement, FastaMap, Nucleotide, Strand};
use std::collections::{HashMap, HashSet};

use crate::read_filtering::soft_clip::is_within_edits_of_a_prefix;
use super::splice_junction_cluster::SpliceJunctionCluster;
use super::junction_placement::{place_read_junctions, pool_insertion_abutting_junctions};
use crate::reference::junction_orientation::{is_aligned_against_transcript, mirror_splice_junction};
use crate::options::{RNAVariantCallingOptions, RNAVariantFilteringOptions};


struct Junction {
    cluster: usize,

    /// The consensus of its records, under the most specific label among them.
    operation: GraphOperation,

    /// Per side, the (model index, record) of the read that spells most past the junction from
    /// it, among the records whose read positions are the junction's flanks; `None` when none is.
    spellers: [Option<(usize, VariantRecord)>; 2]
}


pub(crate) fn resolve_breakpoint_clips(
    models: &mut [TranscriptModel],
    clusters: &mut [SpliceJunctionCluster],
    chromosomes: &BiMap<Box<str>, u16>,
    fasta: &FastaMap,
    calling: &RNAVariantCallingOptions,
    filtering: &RNAVariantFilteringOptions,
    support_floor: impl Fn(u32) -> usize
) -> usize {
    let distance: u32 = calling.soft_clip_max_boundary_distance;
    let mut replacements: Vec<(usize, VariantRecord, VariantRecord)> = Vec::new();
    let mut moves: Vec<(usize, usize, usize)> = Vec::new();
    {
        let models: &[TranscriptModel] = models;
        let index: HashMap<usize, usize> = models
            .iter()
            .enumerate()
            .map(|(i, model)| (model.get_read_id(), i))
            .collect();

        // Step 1. The junctions each provisional cluster supports. A junction arrives under the
        // label each read's annotation gave it, and from reads of either orientation with both
        // strands flipped, so records pool on position and operation alone.
        let mut junctions: Vec<Junction> = Vec::new();
        for (c, cluster) in clusters.iter().enumerate() {
            if cluster.read_ids.len() < filtering.min_total_depth {
                continue;
            }
            let floor: usize = support_floor(cluster.read_ids.len() as u32).max(filtering.min_reads);
            let mut records: Vec<(usize, &VariantRecord)> = cluster
                .read_ids
                .iter()
                .flat_map(|read_id| {
                    let i: usize = index[read_id];
                    models[i].get_variant_records().iter().map(move |record| (i, record))
                })
                .filter(|(_, record)| is_breakpoint(record))
                .collect();
            records.sort_by_key(|(_, record)| (
                record.get_chromosome_1(),
                record.get_position_1(),
                record.get_chromosome_2(),
                record.get_position_2(),
                record.get_read_id()
            ));
            let mut pools: Vec<Vec<(usize, &VariantRecord)>> = Vec::new();
            for record in records {
                match pools.iter_mut().find(|pool| same_junction(pool[0].1.get_graph_operation(), record.1.get_graph_operation(), distance)) {
                    Some(pool) => pool.push(record),
                    None => pools.push(vec![record])
                }
            }
            for pool in pools {
                if pool.iter().map(|(_, record)| record.get_read_id()).collect::<HashSet<usize>>().len() < floor {
                    continue;
                }
                let label: VariantType = pool
                    .iter()
                    .map(|(_, record)| record.get_variant_type())
                    .max_by_key(|variant_type| breakpoint_label_rank(variant_type))
                    .unwrap()
                    .clone();
                let operation: GraphOperation = VariantCall::from_variant_records(
                    0,
                    pool.iter().map(|(_, record)| (*record).clone()).collect(),
                    calling.poa_match_score,
                    calling.poa_mismatch_score,
                    calling.poa_gap_open_score,
                    calling.poa_gap_extend_score
                ).get_consensus_graph_operation().with_variant_type(label);
                // A record's read positions are a junction's flanks when the read aligned across it;
                // a junction retyped from an insertion keeps the insertion's read positions, so its
                // read spells nothing reliable past this side.
                let speller = |s: usize| -> Option<(usize, VariantRecord)> {
                    let reach = |(i, record): &&(usize, &VariantRecord)| -> (u32, std::cmp::Reverse<usize>) {
                        let (flank, step): (u32, i64) = flank(record, s);
                        let length: u32 = if step > 0 { models[*i].get_alignment_model().num_bases() - flank - 1 } else { flank };
                        (length, std::cmp::Reverse(record.get_read_id()))
                    };
                    pool.iter()
                        .filter(|(i, record)| {
                            let (chromosome, position, _, _) = side(record.get_graph_operation(), s);
                            let alignment: &AlignmentModel = models[*i].get_alignment_model();
                            alignment.get_base(flank(record, s).0).get_placement().get_coordinate()
                                .is_some_and(|(c, p, _)| c == chromosome && p == position)
                        })
                        .max_by_key(reach)
                        .map(|(i, record)| (*i, (*record).clone()))
                };
                junctions.push(Junction { cluster: c, operation, spellers: [speller(0), speller(1)] });
            }
        }
        if junctions.is_empty() {
            return 0;
        }

        // Clusters on either side of a fusion each support it: one junction, however many
        // clusters hold it.
        let mut identity: Vec<usize> = Vec::with_capacity(junctions.len());
        for (k, junction) in junctions.iter().enumerate() {
            let first: usize = (0..k)
                .find(|&j| same_junction(&junctions[j].operation, &junction.operation, distance))
                .map_or(k, |j| identity[j]);
            identity.push(first);
        }
        let supports: HashSet<(usize, usize)> = junctions
            .iter()
            .zip(identity.iter())
            .map(|(junction, id)| (junction.cluster, *id))
            .collect();

        // Step 2. Every read's chain as the clusterer keyed it: boundaries placed, transcript
        // orientation.
        let memberships: HashMap<usize, usize> = clusters
            .iter()
            .enumerate()
            .flat_map(|(c, cluster)| cluster.read_ids.iter().map(move |read_id| (*read_id, c)))
            .collect();
        let mut placed = models
            .iter()
            .map(|model| place_read_junctions(model, chromosomes, fasta))
            .collect::<Vec<_>>();
        pool_insertion_abutting_junctions(&mut placed);

        // Step 3. The junctions each read shows, and where it goes.
        for (i, model) in models.iter().enumerate() {
            let Some(&source) = memberships.get(&model.get_read_id()) else {
                continue;
            };
            let mut shown: HashSet<usize> = HashSet::new();
            for record in model.get_variant_records() {
                if is_breakpoint(record) {
                    shown.extend(
                        junctions
                            .iter()
                            .zip(identity.iter())
                            .filter(|(junction, _)| same_junction(&junction.operation, record.get_graph_operation(), distance))
                            .map(|(_, id)| *id)
                    );
                } else if let Some((k, replacement)) = resolve_clip(record, model, models, &junctions, &identity, distance, calling.soft_clip_bases_per_edit, calling.soft_clip_min_partner_bases) {
                    replacements.push((i, record.clone(), replacement));
                    shown.insert(identity[k]);
                }
            }
            if shown.is_empty() {
                continue;
            }
            let mut chain: Vec<SpliceJunction> = placed[i].iter().map(|placed| placed.junction.clone()).collect();
            if is_aligned_against_transcript(&chain, chromosomes, fasta) {
                chain = chain.iter().rev().map(mirror_splice_junction).collect();
            }
            let targets: Vec<usize> = (0..clusters.len())
                .filter(|&c| {
                    shown.iter().all(|id| supports.contains(&(c, *id)))
                        && (chain.is_empty() || clusters[c].splice_junctions.windows(chain.len()).any(|window| window == chain))
                })
                .collect();
            if !targets.contains(&source) && targets.len() == 1 {
                moves.push((model.get_read_id(), source, targets[0]));
            }
        }
    }

    for (i, original, replacement) in replacements.iter() {
        models[*i].replace_variant_record(original, replacement.clone());
    }
    for (read_id, source, target) in moves {
        clusters[source].read_ids.remove(&read_id);
        clusters[target].read_ids.insert(read_id);
    }
    replacements.len()
}


fn is_breakpoint(record: &VariantRecord) -> bool {
    matches!(
        record.get_variant_type(),
        VariantType::Breakpoint | VariantType::FusionGene | VariantType::Translocation
    ) && record.is_resolved()
}


fn same_junction(a: &GraphOperation, b: &GraphOperation, distance: u32) -> bool {
    a.get_chromosome_1() == b.get_chromosome_1()
        && a.get_chromosome_2() == b.get_chromosome_2()
        && a.get_operation_type_1() == b.get_operation_type_1()
        && a.get_operation_type_2() == b.get_operation_type_2()
        && a.get_position_1().abs_diff(b.get_position_1()) <= distance
        && a.get_position_2().abs_diff(b.get_position_2()) <= distance
}


/// The side of an operation: (chromosome, position, operation type, strand).
fn side(operation: &GraphOperation, s: usize) -> (u16, u32, &GraphOperationType, &Strand) {
    if s == 0 {
        (operation.get_chromosome_1(), operation.get_position_1(), operation.get_operation_type_1(), operation.get_strand_1())
    } else {
        (operation.get_chromosome_2(), operation.get_position_2(), operation.get_operation_type_2(), operation.get_strand_2())
    }
}


fn flank(record: &VariantRecord, s: usize) -> (u32, i64) {
    let (_, _, operation, strand) = side(record.get_graph_operation(), s);
    let first: u32 = record.get_read_position_1().min(record.get_read_position_2());
    let second: u32 = record.get_read_position_1().max(record.get_read_position_2());
    if (*operation == GraphOperationType::Downstream) == (*strand == Strand::Forward) {
        (first, 1)
    } else {
        (second, -1)
    }
}


fn spell(alignment: &AlignmentModel, from: u32, step: i64, into: u32, length: usize, complement: bool) -> String {
    let mut position: i64 = from as i64 + step * (1 - into as i64);
    let mut bases: String = String::new();
    while bases.len() < length && 0 <= position && position < alignment.num_bases() as i64 {
        let nucleotide: &Nucleotide = alignment.get_base(position as u32).get_nucleotide();
        bases.push_str(if complement { nucleotide.complement() } else { nucleotide.clone() }.as_str());
        position += step;
    }
    bases
}


fn resolve_clip(
    clip: &VariantRecord,
    model: &TranscriptModel,
    models: &[TranscriptModel],
    junctions: &[Junction],
    identity: &[usize],
    distance: u32,
    bases_per_edit: u32,
    min_partner_bases: u32
) -> Option<(usize, VariantRecord)> {
    if *clip.get_variant_type() != VariantType::Insertion {
        return None;
    }
    let alignment: &AlignmentModel = model.get_alignment_model();
    let start: u32 = clip.get_read_position_1().min(clip.get_read_position_2());
    let end: u32 = clip.get_read_position_1().max(clip.get_read_position_2());
    let leading: bool = start == 0;
    if !leading && end + 1 != alignment.num_bases() {
        return None;
    }
    // The aligned base beside the clip, the step from it into the clip, and the way the read's
    // arm runs from it on the genome.
    let (anchor, away): (u32, i64) = if leading { (end + 1, -1) } else { (start.checked_sub(1)?, 1) };
    if anchor >= alignment.num_bases() {
        return None;
    }
    let (chromosome, position, strand) = alignment.get_base(anchor).get_placement().get_coordinate()?;
    let arm: i64 = if leading == (*strand == Strand::Forward) { 1 } else { -1 };

    let mut resolved: Option<(usize, VariantRecord)> = None;
    let mut ids: HashSet<usize> = HashSet::new();
    for (k, junction) in junctions.iter().enumerate() {
        // The side the read's arm is on: its chromosome, within reach of the anchor, with its arm
        // running the same way (below a Downstream side, above an Upstream one).
        let near: Vec<usize> = (0..2)
            .filter(|&s| {
                let (c, p, operation, _) = side(&junction.operation, s);
                c == chromosome
                    && (p as i64 - position as i64).unsigned_abs() <= distance as u64
                    && (*operation == GraphOperationType::Downstream) == (arm < 0)
            })
            .collect();
        let &[s] = near.as_slice() else {
            continue;
        };

        // The read's clip, and the bases it aligned past the speller's flank, read away from its
        // arm; the speller's read from the clip's first base on.
        let Some((i, speller)) = &junction.spellers[s] else {
            continue;
        };
        let (_, speller_position, _, speller_strand) = side(speller.get_graph_operation(), s);
        let delta: i64 = (position as i64 - speller_position as i64) * arm;
        let query: String = spell(alignment, anchor, away, (-delta).max(0) as u32, usize::MAX, false);
        let edits: usize = if bases_per_edit == 0 { 0 } else { query.len() / bases_per_edit as usize };
        let (from, step): (u32, i64) = flank(speller, s);
        let target: String = spell(models[*i].get_alignment_model(), from, step, delta.max(0) as u32, query.len() + edits, speller_strand != strand);
        let untemplated: usize = speller.get_read_position_1().abs_diff(speller.get_read_position_2()).saturating_sub(1) as usize;
        if query.len().saturating_sub(delta.max(0) as usize + untemplated) < min_partner_bases as usize
            || !is_within_edits_of_a_prefix(query.as_bytes(), target.as_bytes(), query.len() - edits, (query.len() + edits).min(target.len()), edits) {
            continue;
        }
        ids.insert(identity[k]);
        if resolved.is_none() {
            resolved = junction_record(clip, alignment.num_bases(), anchor, away, arm, position, strand, &junction.operation, s).map(|record| (k, record));
        }
    }
    if ids.len() == 1 { resolved } else { None }
}


fn junction_record(
    clip: &VariantRecord,
    read_length: u32,
    anchor: u32,
    away: i64,
    arm: i64,
    position: u32,
    strand: &Strand,
    operation: &GraphOperation,
    s: usize
) -> Option<VariantRecord> {
    let (_, junction_position, _, junction_strand) = side(operation, s);
    let delta: i64 = (position as i64 - junction_position as i64) * arm;
    let near: i64 = anchor as i64 + away * delta;
    let far: i64 = near + away * (operation.get_sequence_length() as i64 + 1);
    if near.min(far) < 0 || near.max(far) >= read_length as i64 {
        return None;
    }
    let flip: bool = junction_strand != strand;
    let orient = |strand: &Strand| -> Strand {
        match (flip, strand) {
            (true, Strand::Forward) => Strand::Reverse,
            (true, Strand::Reverse) => Strand::Forward,
            _ => strand.clone()
        }
    };
    Some(VariantRecord::new(
        clip.get_read_id(),
        near.min(far) as u32,
        near.max(far) as u32,
        GraphOperation::new(
            operation.get_chromosome_1(),
            operation.get_position_1(),
            orient(operation.get_strand_1()),
            operation.get_operation_type_1().clone(),
            operation.get_chromosome_2(),
            operation.get_position_2(),
            orient(operation.get_strand_2()),
            operation.get_operation_type_2().clone(),
            if flip { reverse_complement(operation.get_sequence()) } else { operation.get_sequence().into() },
            operation.get_variant_type().clone()
        )
    ))
}


#[cfg(test)]
#[path = "../tests/clustering/breakpoint_clip.rs"]
mod tests;
