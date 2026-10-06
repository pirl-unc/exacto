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
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::clustering::splice_junction_cluster::SpliceJunctionCluster;
use crate::read_cluster::rna_read_cluster::RNAReadCluster;
use crate::read_cluster::rna_read_cluster_set::{FailedRNAVariantCall, RNAReadClusterSet};
use crate::variant_filtering::cluster_depths::ClusterDepths;
use crate::variant_filtering::variant_judge::min_reads_at_site;


pub(crate) fn store_cells(
    splice_junction_cluster: &SpliceJunctionCluster,
    variant_calls: &[VariantCall],
    genotypes: HashMap<usize, Vec<(usize, Allele)>>,
    failed_variant_calls: Vec<FailedVariantCall>,
    (cell_read_ids, shared_read_ids): (HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>),
    repeat_lengths: &HashMap<usize, Option<u32>>,
    transcript_models_map: &HashMap<usize, Arc<TranscriptModel>>,
    read_support_index: &RNAVariantReadSupportIndex,
    next_cluster_id: &mut usize,
    rna_read_cluster_set: &mut RNAReadClusterSet
) {
    let mut phased_cluster_id: usize = *next_cluster_id;
    let reference_gene_transcript_ids: HashSet<(Box<str>, Box<str>)> = splice_junction_cluster
        .reference_gene_transcript_ids
        .iter()
        .cloned()
        .collect();

    // The output cluster ID of every cell, in cell order, and the output cluster every
    // read was assigned to. A read shared with other cells is filed under its own.
    let mut cell_ids: Vec<usize> = cell_read_ids.keys().copied().collect();
    cell_ids.sort_unstable();
    // HashMap<read ID, output cluster ID>
    let output_cluster_of_read: HashMap<usize, usize> = cell_ids
        .iter()
        .enumerate()
        .flat_map(|(index, cell_id)| {
            cell_read_ids[cell_id].iter().map(move |read_id| (*read_id, phased_cluster_id + index))
        })
        .collect();
    // The output cluster holding most of `read_ids`, the lowest among equals, and how many.
    let holding_most = |read_ids: &[usize]| -> (usize, u32) {
        let mut counts: BTreeMap<usize, u32> = BTreeMap::new();
        for read_id in read_ids.iter() {
            *counts.entry(output_cluster_of_read[read_id]).or_insert(0) += 1;
        }
        counts
            .into_iter()
            .max_by(|(id_a, count_a), (id_b, count_b)| count_a.cmp(count_b).then_with(|| id_b.cmp(id_a)))
            .unwrap()
    };
    
    let members = |cell_id: &usize| -> HashSet<usize> {
        let mut read_ids: HashSet<usize> = cell_read_ids[cell_id].clone();
        read_ids.extend(shared_read_ids.get(cell_id).into_iter().flatten().copied());
        read_ids
    };
    
    // Vec<(output cluster ID, reads of the cell, depths of the cell)>
    let cells: Vec<(usize, HashSet<usize>, ClusterDepths)> = cell_ids
        .iter()
        .enumerate()
        .map(|(index, cell_id)| {
            let read_ids: HashSet<usize> = members(cell_id);
            let transcript_models: Vec<Arc<TranscriptModel>> = read_ids
                .iter()
                .map(|read_id| Arc::clone(&transcript_models_map[read_id]))
                .collect();
            (phased_cluster_id + index, read_ids, ClusterDepths::new(&transcript_models))
        })
        .collect();
    // HashMap<output cluster ID, HashSet<variant call ID>>
    let mut carried: HashMap<usize, HashSet<usize>> = HashMap::new();
    for variant_call in variant_calls.iter() {
        let op: &GraphOperation = variant_call.get_consensus_graph_operation();
        let call_read_ids: Vec<usize> = variant_call.get_read_ids();
        let repeat_length: Option<u32> = repeat_lengths[&variant_call.get_id()];
        // (reads of the call in the cell, the cell's depth at the site, the floor there)
        let tally = |read_ids: &HashSet<usize>, depths: &ClusterDepths| -> (u32, u32, usize) {
            let carriers: Vec<&TranscriptModel> = call_read_ids
                .iter()
                .filter(|read_id| read_ids.contains(read_id))
                .map(|read_id| &*transcript_models_map[read_id])
                .collect();
            let num_reads: u32 = carriers.len() as u32;
            let (depth_1, depth_2): (u32, u32) = depths.at_site(op, &carriers);
            let min_reads: usize = repeat_length.map_or(1, |repeat_length| {
                min_reads_at_site(read_support_index, repeat_length, (depth_1, depth_2)) as usize
            });
            (num_reads, depth_1.max(depth_2), min_reads)
        };
        let mut is_carried: bool = false;
        for (output_cluster_id, read_ids, depths) in cells.iter() {
            let (num_reads, _, min_reads): (u32, u32, usize) = tally(read_ids, depths);
            if num_reads > 0 && num_reads as usize >= min_reads {
                carried.entry(*output_cluster_id).or_default().insert(variant_call.get_id());
                is_carried = true;
            }
        }
        if !is_carried {
            let (output_cluster_id, _): (usize, u32) = holding_most(&call_read_ids);
            let (_, read_ids, depths) = cells.iter().find(|(id, _, _)| *id == output_cluster_id).unwrap();
            let (num_reads, depth, min_reads): (u32, u32, usize) = tally(read_ids, depths);
            rna_read_cluster_set.add_failed_variant_call(FailedRNAVariantCall {
                cluster_id: output_cluster_id,
                graph_operation: op.clone(),
                num_reads,
                total_depth: depth as i32,
                failure: VariantCallFailure::TooFewReads {
                    num_reads,
                    min_reads: u32::try_from(min_reads).unwrap_or(u32::MAX),
                    repeat_len: repeat_length
                }
            });
        }
    }
    drop(cells);

    // The calls the gates rejected, each given to the output cluster holding most of its
    // reads.
    for failed in failed_variant_calls {
        let (output_cluster_id, _): (usize, u32) = holding_most(&failed.variant_call.get_read_ids());
        rna_read_cluster_set.add_failed_variant_call(FailedRNAVariantCall {
            cluster_id: output_cluster_id,
            graph_operation: failed.variant_call.get_consensus_graph_operation().clone(),
            num_reads: failed.variant_call.get_num_reads() as u32,
            total_depth: failed.variant_call.get_total_depth(),
            failure: failed.failure
        });
    }

    for cell_id in cell_ids.iter() {
        let read_ids: HashSet<usize> = members(cell_id);
        let carried_ids: HashSet<usize> = carried.remove(&phased_cluster_id).unwrap_or_default();

        // The calls this cell carries, and the cell's own rows of each.
        let cell_variant_calls: Vec<VariantCall> = variant_calls
            .iter()
            .filter(|variant_call| carried_ids.contains(&variant_call.get_id()))
            .cloned()
            .collect();
        // HashMap<variant call ID, HashMap<read ID, Allele>>
        let cell_genotypes: HashMap<usize, HashMap<usize, Allele>> = genotypes
            .iter()
            .filter(|(variant_call_id, _)| carried_ids.contains(variant_call_id))
            .map(|(variant_call_id, alleles)| {
                let rows: HashMap<usize, Allele> = alleles
                    .iter()
                    .filter(|(read_id, _)| read_ids.contains(read_id))
                    .map(|(read_id, allele)| (*read_id, allele.clone()))
                    .collect();
                (*variant_call_id, rows)
            })
            .collect();

        rna_read_cluster_set.add_cluster(RNAReadCluster::new(
            phased_cluster_id,
            read_ids,
            splice_junction_cluster.splice_junctions.clone(),
            cell_variant_calls,
            cell_genotypes,
            reference_gene_transcript_ids.clone()
        ).with_shared_read_ids(shared_read_ids.get(cell_id).cloned().unwrap_or_default()));
        phased_cluster_id += 1;
    }
    *next_cluster_id = phased_cluster_id;
}
