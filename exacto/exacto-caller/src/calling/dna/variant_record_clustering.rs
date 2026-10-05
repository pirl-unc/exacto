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
use edit_distance::edit_distance;
use exacto_core::prelude::*;
use rayon::prelude::*;
use rayon::iter::IntoParallelIterator;
use std::collections::{HashMap, HashSet};
use std::ops::RangeInclusive;
use std::sync::Arc;

use crate::prelude::*;
use crate::calling::dna::breakpoint_rescue::retype_insertion_as_breakpoints;


pub(crate) fn build_dna_breakpoint_interval_trees(
    variant_records: &Vec<Arc<VariantRecord>>,
    breakend_rescue_search_distance: u32
) -> HashMap<ReferenceChromosomeID, IntervalTree<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)>> {
    // HashMap<chromosome, <IntervalTree<(chromosome 1, position 1, chromosome 2, position 2, variant record index)>>
    let mut interval_tree_map: HashMap<ReferenceChromosomeID, IntervalTree<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)>> = HashMap::new();
    for (i, variant_record) in variant_records.iter().enumerate() {
        if variant_record.get_variant_type() == &VariantType::Breakpoint
            || variant_record.get_variant_type() == &VariantType::Translocation {
            let interval_1: Interval<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)> = Interval::new(
                variant_record.get_position_1() as isize - breakend_rescue_search_distance as isize,
                variant_record.get_position_1() as isize + breakend_rescue_search_distance as isize,
                (
                    variant_record.get_chromosome_1(),
                    variant_record.get_position_1(),
                    variant_record.get_chromosome_2(),
                    variant_record.get_position_2(),
                    i
                )
            );
            let interval_2: Interval<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)> = Interval::new(
                variant_record.get_position_2() as isize - breakend_rescue_search_distance as isize,
                variant_record.get_position_2() as isize + breakend_rescue_search_distance as isize,
                (
                    variant_record.get_chromosome_1(),
                    variant_record.get_position_1(),
                    variant_record.get_chromosome_2(),
                    variant_record.get_position_2(),
                    i
                )
            );
            interval_tree_map
                .entry(variant_record.get_chromosome_1())
                .or_insert(IntervalTree::new())
                .insert(interval_1);
            interval_tree_map
                .entry(variant_record.get_chromosome_2())
                .or_insert(IntervalTree::new())
                .insert(interval_2);
        }
    }
    interval_tree_map
}


pub(crate) fn calculate_max_dna_clustering_distance(
    size: u32,
    sequencing_error: f64,
    max_distance: u32
) -> u32 {
    assert!(
        (0.0..=1.0).contains(&sequencing_error),
        "sequencing_error must be in [0.0, 1.0], got {sequencing_error}"
    );
    if max_distance == 0 {
        return 0;
    }
    let max_distance_: f64 = max_distance as f64;
    let numerator: f64 = -sequencing_error * size as f64 * max_distance_.sqrt();
    (max_distance as f64 * (1.0 - f64::exp(numerator / max_distance_))).ceil() as u32
}


pub(crate) fn cluster_dna_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    owned_positions: RangeInclusive<ReferencePosition>,
    fasta_map: &FastaMap,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    num_threads: usize,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64,
    bkpt_rescue: bool,
    bkpt_rescue_min_query_length: u32,
    bkpt_rescue_max_query_length: u32,
    bkpt_rescue_search_distance: u32,
    bkpt_rescue_gap_open: i32,
    bkpt_rescue_gap_extend: i32,
    bkpt_rescue_k: u32,
    bkpt_rescue_band_width: u32,
    bkpt_rescue_min_score_fraction: f64,
    bkpt_rescue_min_query_coverage: f64,
    bkpt_rescue_min_span_proportion: f64,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32
) -> Vec<VariantCall> {
    // Step 1. Create an IntervalTree of breakpoint (BND and TRA) positions.
    // HashMap<reference chromosome ID, IntervalTree<(reference chromosome ID, reference position)>>
    let mut interval_tree_map: HashMap<ReferenceChromosomeID, IntervalTree<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)>> = build_dna_breakpoint_interval_trees(
        &variant_records,
        bkpt_rescue_search_distance
    );

    // Step 2. Retype insertion sequences as breakpoints.
    // Try to align insertion sequences to the BND and TRA mate positions.
    // HashMap<retyped variant record, position 1 of the insertion it was retyped from>
    let mut retyped_positions: HashMap<VariantRecord, ReferencePosition> = HashMap::new();
    let variant_records: Vec<Arc<VariantRecord>> = variant_records
        .into_iter()
        .flat_map(|vr| {
            // Only retype insertions, and only when rescue is on.
            if bkpt_rescue == false || vr.get_variant_type() != &VariantType::Insertion {
                return vec![vr];
            }

            // If the interval trees are empty for the insertion chromosome, then don't realign.
            let Some(interval_tree) = interval_tree_map.get(&vr.get_chromosome_1()) else {
                return vec![vr];
            };

            // Get mate 2 positions of BND or TRA whose mate 1 position is near the insertion position.
            let breakpoint_positions: Vec<&(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)> = interval_tree
                .overlaps(vr.get_position_1() as isize, vr.get_position_2() as isize);

            if breakpoint_positions.is_empty() {
                // No BND or TRA exists near the insertion.
                vec![vr]
            } else {
                // Align the insertion sequence.
                let mut query_positions: HashMap<ReferenceChromosomeID, Vec<ReferencePosition>> = HashMap::new();
                for (chromosome_1, position_1, chromosome_2, position_2, variant_record_index) in breakpoint_positions.iter() {
                    query_positions
                        .entry(*chromosome_1)
                        .or_insert_with(Vec::new)
                        .push(*position_1);
                    query_positions
                        .entry(*chromosome_2)
                        .or_insert_with(Vec::new)
                        .push(*position_2);
                }
                match retype_insertion_as_breakpoints(
                    vr.as_ref(),
                    &query_positions,
                    fasta_map,
                    chromosome_names_map,
                    bkpt_rescue_min_query_length,
                    bkpt_rescue_max_query_length,
                    bkpt_rescue_search_distance,
                    bkpt_rescue_gap_open,
                    bkpt_rescue_gap_extend,
                    bkpt_rescue_k,
                    bkpt_rescue_band_width,
                    bkpt_rescue_min_score_fraction,
                    bkpt_rescue_min_query_coverage,
                    bkpt_rescue_min_span_proportion
                ) {
                    Some(retyped_records) => {
                        for retyped_record in retyped_records.iter() {
                            retyped_positions.insert(retyped_record.clone(), vr.get_position_1());
                        }
                        retyped_records
                            .into_iter()
                            .map(Arc::new)
                            .collect::<Vec<_>>()
                    },
                    None => {
                        vec![vr]
                    }
                }
            }
        })
        .collect();

    // Step 3. Split the variant records into breakpoint records and non-breakpoint records.
    let (breakpoint_records, nonbreakpoint_records): (Vec<Arc<VariantRecord>>, Vec<Arc<VariantRecord>>) =
        variant_records
            .into_iter()
            .map(|rc| Arc::new((*rc).clone()))
            .partition(|vr| {
                matches!(
                    vr.get_variant_type(),
                    VariantType::Breakpoint | VariantType::Translocation
                )
            });

    // Step 4. Cluster breakpoint variant records.
    let mut clusters: Vec<VariantRecordCluster> = cluster_breakpoint_dna_variant_records(
        breakpoint_records,
        max_clustering_distance,
        min_size_proportion,
        max_ins_norm_edit_distance,
        sequencing_error,
        num_threads
    );

    // Step 5. Cluster non-breakpoint (local) variant records.
    let variant_records_map: HashMap<(ReferenceChromosomeID, ReferenceChromosomeID, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> = split_variant_records(
        nonbreakpoint_records,
        true,
        num_threads
    );
    for curr_variant_records in variant_records_map.values() {
        clusters.extend(cluster_non_breakpoint_dna_variant_records(
            curr_variant_records,
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_clustering_distance,
            sequencing_error,
            num_threads
        ));
    }

    // Step 6. Create the variant calls of the clusters that start in `owned_positions`.
    clusters.retain(|cluster| {
        let start: ReferencePosition = cluster
            .get_variant_records()
            .iter()
            .map(|variant_record| *retyped_positions.get(variant_record.as_ref()).unwrap_or(&variant_record.get_position_1()))
            .min()
            .unwrap();
        owned_positions.contains(&start)
    });
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let variant_calls: Vec<VariantCall> = thread_pool.install(|| {
        clusters
            .par_chunks(((clusters.len() + num_threads - 1) / num_threads).max(1))
            .flat_map(|variant_record_clusters| {
                let mut variant_calls: Vec<VariantCall> = Vec::new();
                for variant_record_cluster in variant_record_clusters.iter() {
                    // Get the variant records from the cluster.
                    let variant_records: HashSet<VariantRecord> = variant_record_cluster
                        .get_variant_records()
                        .iter()
                        .map(Arc::as_ref)
                        .cloned()
                        .collect();

                    // Create a variant call.
                    let variant_call: VariantCall = VariantCall::from_variant_records(
                        0,
                        variant_records,
                        match_score,
                        mismatch_score,
                        gap_open_score,
                        gap_extend_score,
                    );

                    variant_calls.push(variant_call);
                }
                variant_calls.into_par_iter()
            })
            .collect()
    });

    variant_calls
}


fn cluster_breakpoint_dna_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    max_clustering_distance: u32,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    sequencing_error: f64,
    num_threads: usize
) -> Vec<VariantRecordCluster> {
    assert!(
        variant_records.iter().all(|variant_record| matches!(
            variant_record.get_variant_type(),
            VariantType::Breakpoint | VariantType::Translocation
        ))
    );

    if variant_records.is_empty() {
        return vec![];
    }

    // Check if a variant record has resolved graph operations.
    let is_resolved = |i: usize| -> bool {
        *variant_records[i].get_operation_2() != GraphOperationType::Noop
    };

    // Step 1. Group variant records by chromosome 1 and operation 1.
    let mut position_1_groups_map: HashMap<(ReferenceChromosomeID, GraphOperationType), Vec<usize>> = HashMap::new();
    for (i, variant_record) in variant_records.iter().enumerate() {
        position_1_groups_map
            .entry((
                variant_record.get_chromosome_1(),
                variant_record.get_operation_1().clone(),
            ))
            .or_default()
            .push(i);
    }

    // Step 2. Group variant records by chromosome 2 and operation 2.
    let mut position_2_groups_map: HashMap<(ReferenceChromosomeID, GraphOperationType), Vec<usize>> = HashMap::new();
    for (i, variant_record) in variant_records.iter().enumerate() {
        if is_resolved(i) {
            position_2_groups_map
                .entry((
                    variant_record.get_chromosome_2(),
                    variant_record.get_operation_2().clone(),
                ))
                .or_default()
                .push(i);
        }
    }

    // Step 3. Sort each group by position 1.
    for group in position_1_groups_map.values_mut() {
        group.sort_unstable_by_key(|&i| variant_records[i].get_position_1());
    }

    // Step 4. Sort each group by position 2.
    for group in position_2_groups_map.values_mut() {
        group.sort_unstable_by_key(|&i| variant_records[i].get_position_2());
    }

    // Step 5. Start a thread pool.
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    // Step 6. Cluster variant records by chromosome 1, operation 1.
    let merge_by_position_1 = |indices: &[usize]| -> Vec<HashSet<usize>> {
        let mut uf: UnionFind = UnionFind::new();
        for (k, &i) in indices.iter().enumerate() {
            uf.union(i as u32, i as u32);
            let a: &VariantRecord = variant_records[i].as_ref();
            for &j in &indices[k + 1..] {
                let b: &VariantRecord = variant_records[j].as_ref();
                if b.get_position_1() - a.get_position_1() > max_clustering_distance {
                    break;
                }
                uf.union(i as u32, j as u32);
            }
        }
        uf.get_clusters()
            .into_iter()
            .map(|cluster| cluster.into_iter().map(|i| i as usize).collect())
            .collect()
    };
    let position_1_groups_list: Vec<HashSet<usize>> = thread_pool.install(|| {
        position_1_groups_map
            .par_iter()
            .flat_map_iter(|(_, indices)| merge_by_position_1(indices))
            .collect()
    });

    // Step 7. Cluster variant records by chromosome 2, operation 2.
    let merge_by_position_2 = |indices: &[usize]| -> Vec<HashSet<usize>> {
        let mut uf: UnionFind = UnionFind::new();
        for (k, &i) in indices.iter().enumerate() {
            uf.union(i as u32, i as u32);
            let a: &VariantRecord = variant_records[i].as_ref();
            for &j in &indices[k + 1..] {
                let b: &VariantRecord = variant_records[j].as_ref();
                if b.get_position_2() - a.get_position_2() > max_clustering_distance {
                    break;
                }
                uf.union(i as u32, j as u32);
            }
        }
        uf.get_clusters()
            .into_iter()
            .map(|cluster| cluster.into_iter().map(|i| i as usize).collect())
            .collect()
    };
    let position_2_groups_list: Vec<HashSet<usize>> = thread_pool.install(|| {
        position_2_groups_map
            .par_iter()
            .flat_map_iter(|(_, indices)| merge_by_position_2(indices))
            .collect()
    });

    // Step 8. Assign an ID to groups 1 and 2.
    let mut group_id: usize = 1; // reserve 0 for orphaned insertion variant records' group 2 ID
    // HashMap<variant record index, group 1 ID>
    let mut group_1_ids_map: HashMap<usize, usize> = HashMap::new();
    // HashMap<variant record index, group 1 ID>
    let mut group_2_ids_map: HashMap<usize, usize> = HashMap::new();
    for indices in position_1_groups_list.iter() {
        for &i in indices.iter() {
            group_1_ids_map.insert(i, group_id);
        }
        group_id += 1;
    }
    for indices in position_2_groups_list.iter() {
        for &i in indices.iter() {
            group_2_ids_map.insert(i, group_id);
        }
        group_id += 1;
    }

    // Step 9. Merge resolved variant records based on their groups 1 and 2 IDs.
    //
    // If two variant records (A and B) are resolved, they can be clustered if BOTH breakpoints agree:
    //
    // |position_1(A) - position_1(B)| <= max_clustering_distance AND
    // |position_2(A) - position_2(B)| <= max_clustering_distance
    //
    // Additionally, A and B must have the same
    // (chromosome_1, chromosome_2) and (operation_1, operation_2).

    // HashMap<(group 1 ID, group 2 ID), Vec<variant record index>>
    let mut groups: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for i in 0..variant_records.len() {
        if is_resolved(i) {
            let group_1_id: usize = *group_1_ids_map.get(&i).unwrap();
            let group_2_id: usize = *group_2_ids_map.get(&i).unwrap();
            groups.entry((group_1_id, group_2_id)).or_insert_with(Vec::new).push(i);
        }
    }

    // Step 10. Create a map of group 1 ID to group 2 IDs.
    // HashMap<group 1 ID, HashSet<group 2 IDs>>
    let mut group_1_group_2_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    for (group_1_id, group_2_id) in groups.keys() {
        group_1_group_2_map.entry(*group_1_id).or_insert_with(HashSet::new).insert(*group_2_id);
    }

    // Step 11. Merge unresolved variant records based on their group 1.
    //
    // If one variant record is resolved and the other is unresolved, they can be clustered if
    // position 1 breakpoints agree:
    //
    // |position_1(A) - position_1(B)| <= max_clustering_distance
    let side_2_group = |i: usize| -> Option<(usize, usize)> {
        let clip: &VariantRecord = variant_records[i].as_ref();
        let indices: &Vec<usize> = position_2_groups_map.get(&(clip.get_chromosome_1(), clip.get_operation_1().clone()))?;
        let position: ReferencePosition = clip.get_position_1();
        let max_distance: u32 = calculate_max_dna_clustering_distance(
            clip.get_graph_operation().get_sequence_length() as u32,
            sequencing_error,
            max_clustering_distance
        );
        let from: usize = indices.partition_point(|&j| variant_records[j].get_position_2() < position.saturating_sub(max_distance));
        let to: usize = indices.partition_point(|&j| variant_records[j].get_position_2() <= position.saturating_add(max_distance));
        let junctions: HashSet<(usize, usize)> = indices[from..to]
            .iter()
            .map(|j| (group_1_ids_map[j], group_2_ids_map[j]))
            .collect();
        if junctions.len() == 1 {
            junctions.into_iter().next()
        } else {
            None
        }
    };
    for i in 0..variant_records.len() {
        if !is_resolved(i) {
            let group_1_id: usize = *group_1_ids_map.get(&i).unwrap();
            if group_1_group_2_map.contains_key(&group_1_id) {
                let group_2_ids: &HashSet<usize> = group_1_group_2_map.get(&group_1_id).unwrap();
                if group_2_ids.len() == 1 {
                    let group_2_id: usize = group_2_ids.iter().next().unwrap().clone();
                    groups.entry((group_1_id, group_2_id)).or_insert_with(Vec::new).push(i);
                } else {
                    // Multiple group 2 IDs exist.
                    // Do not cluster this unresolved variant record to any cluster.
                }
            } else if let Some(junction) = side_2_group(i) {
                groups.get_mut(&junction).unwrap().push(i);
            } else {
                // Noop insertion sequence.
                // Set group 2 ID as 0.
                groups.entry((group_1_id, 0)).or_insert_with(Vec::new).push(i);
            }
        }
    }

    // Step 12. Create variant record clusters.
    let mut clusters: Vec<VariantRecordCluster> = Vec::new();
    for ((_, group_2_id), indices) in groups.iter() {
        if *group_2_id == 0 && indices.len() > 1 {
            // Orphan pool: no member resolved a mate. If the clipped tails agree on one
            // inserted allele this is a terminal insertion, not a junction.
            let as_insertions: Vec<Arc<VariantRecord>> =
                rebuild_orphan_pool_as_insertions(&variant_records, indices);
            let insertion_clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_dna_variant_records(
                &as_insertions,
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_clustering_distance,
                sequencing_error,
                num_threads
            );
            if insertion_clusters.len() == 1 {
                clusters.extend(insertion_clusters);
                continue;
            }
        }
        let variant_records_selected: Vec<Arc<VariantRecord>> = indices
            .iter()
            .map(|&i| Arc::clone(&variant_records[i]))
            .collect();
        clusters.push(VariantRecordCluster::from_variant_records(&variant_records_selected));
    }

    clusters
}


fn cluster_non_breakpoint_dna_variant_records(
    variant_records: &Vec<Arc<VariantRecord>>,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64,
    num_threads: usize
) -> Vec<VariantRecordCluster> {
    assert!(
        variant_records.iter().all(|variant_record| !matches!(
            variant_record.get_variant_type(),
            VariantType::Breakpoint | VariantType::Translocation
        ))
    );

    let n: usize = variant_records.len();

    if n == 0 {
        return Vec::new();
    }

    // Step 1. Ensure all records are on the same chromosomes.
    let chromosome_1: ReferenceChromosomeID = variant_records[0].get_chromosome_1();
    let chromosome_2: ReferenceChromosomeID = variant_records[0].get_chromosome_2();
    for vr in variant_records.iter() {
        assert!(
            vr.get_chromosome_1() == chromosome_1,
            "Supplied variant_records must have the same chromosome_1 value."
        );
        assert!(
            vr.get_chromosome_2() == chromosome_2,
            "Supplied variant_records must have the same chromosome_2 value."
        );
    }

    // Step 2. Sort the variant records by position 1, ascending.
    let variant_records: Vec<Arc<VariantRecord>> = {
        let mut sorted: Vec<Arc<VariantRecord>> = variant_records.to_vec();
        sorted.sort_by_key(|variant_record| variant_record.get_position_1());
        sorted
    };

    // Step 3. Identify clusterable pairs with a sliding window over sorted positions.
    let n: usize = variant_records.len();
    let indices: Vec<usize> = (0..variant_records.len()).collect();
    let mates_of = |i: usize| -> (u32, Vec<u32>) {
        let mut mate_ids: Vec<u32> = Vec::new();
        let variant_record_i: &Arc<VariantRecord> = &variant_records[i];
        for j in (i + 1)..n {
            let variant_record_j: &Arc<VariantRecord> = &variant_records[j];
            if meet_dna_clustering_criteria(
                variant_record_i.as_ref(),
                variant_record_j.as_ref(),
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_clustering_distance,
                sequencing_error
            ) {
                mate_ids.push(j as u32);
            }
            if variant_record_j
                .get_position_1()
                .saturating_sub(variant_record_i.get_position_1()) > max_clustering_distance {
                break;
            }
        }
        (i as u32, mate_ids)
    };
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let pairs_list: Vec<(u32, Vec<u32>)> = thread_pool.install(||
        indices
            .par_iter()
            .map(|&i| mates_of(i))
            .collect()
    );

    // Step 4. Union-Find over the mate indices.
    let mut uf: UnionFind = UnionFind::new();
    for (i, mate_indices) in pairs_list.iter() {
        uf.union(*i, *i);
        for j in mate_indices.iter() {
            uf.union(*i, *j);
        }
    }

    // Step 5. Build clusters from the Union-Find sets.
    let mut clusters: Vec<VariantRecordCluster> = Vec::new();
    for indices in uf.get_clusters().iter() {
        let variant_records_selected: Vec<Arc<VariantRecord>> = indices
            .iter()
            .map(|&i| Arc::clone(&variant_records[i as usize]))
            .collect();
        let cluster: VariantRecordCluster = VariantRecordCluster::from_variant_records(
            &variant_records_selected
        );
        clusters.push(cluster);
    }

    clusters
}


pub(crate) fn count_dna_variant_depths(
    variant_calls: &mut [VariantCall],
    bam_file: &str,
    bai_file: &str,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    max_merge_distance: u32
) -> BAMReadDepths {
    let mut positions: HashMap<ReferenceChromosomeName, Vec<ReferencePosition>> = HashMap::new();
    for variant_call in variant_calls.iter() {
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        let chromosome_1: &ReferenceChromosomeName = chromosome_names_map.get_by_right(&operation.get_chromosome_1()).unwrap();
        let chromosome_2: &ReferenceChromosomeName = chromosome_names_map.get_by_right(&operation.get_chromosome_2()).unwrap();
        let (position_1, position_2): (ReferencePosition, ReferencePosition) = (operation.get_position_1(), operation.get_position_2());
        positions.entry(chromosome_1.clone()).or_default().push(position_1);
        positions.entry(chromosome_2.clone()).or_default().push(position_2);
        if chromosome_1 == chromosome_2
            && position_2 > position_1
            && position_2 - position_1 - 1 <= operation.get_sequence_length() as u32 {
            positions.entry(chromosome_1.clone()).or_default().extend(position_1 + 1..position_2);
        }
    }

    // Positions less than `max_merge_distance` bases apart are counted from one query of the BAM
    // file. Only the number of queries depends on it, not the counts.
    let read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, bai_file, &positions, max_merge_distance);

    for variant_call in variant_calls.iter_mut() {
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        let total_depth: ReadDepth = get_dna_variant_position_total_depth(
            &read_depths,
            chromosome_names_map.get_by_right(&operation.get_chromosome_1()).unwrap(),
            operation.get_position_1(),
            operation.get_operation_type_1(),
            chromosome_names_map.get_by_right(&operation.get_chromosome_2()).unwrap(),
            operation.get_position_2(),
            operation.get_operation_type_2(),
            operation.get_standardized_sequence().as_str()
        );
        variant_call.set_total_depth(total_depth as i32);
    }

    read_depths
}


fn get_dna_variant_position_total_depth(
    read_depths: &BAMReadDepths,
    chromosome_1: &str,
    position_1: ReferencePosition,
    operation_1: &GraphOperationType,
    chromosome_2: &str,
    position_2: ReferencePosition,
    operation_2: &GraphOperationType,
    sequence: &str
) -> ReadDepth {
    // SNV
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) == 1 &&
        sequence.len() == 1 {
        return read_depths.get_depth(chromosome_1, position_1 + 1);
    }

    // MNV
    if chromosome_1 == chromosome_2 &&
        *operation_1 == GraphOperationType::Downstream &&
        *operation_2 == GraphOperationType::Upstream &&
        (position_2 - position_1 - 1) == sequence.len() as u32 &&
        sequence.len() >= 2 {
        return (position_1 + 1..position_2)
            .map(|position| read_depths.get_depth(chromosome_1, position))
            .max()
            .unwrap();
    }

    // INS, DEL and BND: the deeper flank.
    read_depths.get_depth(chromosome_1, position_1).max(read_depths.get_depth(chromosome_2, position_2))
}


fn meet_dna_clustering_criteria(
    a: &VariantRecord,
    b: &VariantRecord,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64
) -> bool {
    // Records from the same read ID cannot be clustered.
    let a_read_id: ReadID = a.get_read_id();
    let b_read_id: ReadID = b.get_read_id();
    if a_read_id == b_read_id {
        return false;
    }

    // Variant types must be the same.
    let a_variant_type: &VariantType = a.get_variant_type();
    let b_variant_type: &VariantType = b.get_variant_type();
    if a_variant_type != b_variant_type {
        return false;
    }

    // Orientations must be the same.
    let a_op1: &GraphOperationType = a.get_operation_1();
    let a_op2: &GraphOperationType = a.get_operation_2();
    let b_op1: &GraphOperationType = b.get_operation_1();
    let b_op2: &GraphOperationType = b.get_operation_2();
    if a_op1 != b_op1 || a_op2 != b_op2 {
        return false;
    }

    // Chromosomes 1 and 2 must be the same.
    let a_chr1: ReferenceChromosomeID = a.get_chromosome_1();
    let a_chr2: ReferenceChromosomeID = a.get_chromosome_2();
    let b_chr1: ReferenceChromosomeID = b.get_chromosome_1();
    let b_chr2: ReferenceChromosomeID = b.get_chromosome_2();
    if a_chr1 != b_chr1 || a_chr2 != b_chr2 {
        return false;
    }

    // Strand pairing must be the same.
    let a_s1: &Strand = a.get_strand_1();
    let a_s2: &Strand = a.get_strand_2();
    let b_s1: &Strand = b.get_strand_1();
    let b_s2: &Strand = b.get_strand_2();
    match (a_s1, a_s2) {
        (Strand::Forward, Strand::Forward) => {
            if b_s1 != b_s2 {
                return false;
            }
        },
        (Strand::Forward, Strand::Reverse) => {
            if b_s1 == b_s2 {
                return false;
            }
        },
        (Strand::Reverse, Strand::Forward) => {
            if b_s1 == b_s2 {
                return false;
            }
        },
        (Strand::Reverse, Strand::Reverse) => {
            if b_s1 != b_s2 {
                return false;
            }
        },
        _ => {
            panic!("Invalid strand combination: {} and {}", a.get_strand_1().as_str(), a.get_strand_2().as_str());
        }
    }

    // The size proportion must be within the allowed limit.
    let a_size: isize = a.get_variant_size();
    let b_size: isize = b.get_variant_size();
    let has_valid_sizes: bool = a_size >= 0 && b_size >= 0;
    if has_valid_sizes {
        let min_size: f64 = a_size.min(b_size) as f64;
        let max_size: f64 = a_size.max(b_size) as f64;
        if max_size > 0.0 {
            let size_proportion: f64 = min_size / max_size;
            if size_proportion < min_size_proportion {
                return false;
            }
        }
    }

    // The breakpoint distances must be close where proximity is a function of the variant size.
    let max_distance: u32 = if !has_valid_sizes {
        max_clustering_distance
    } else {
        if a_chr1 == b_chr2 {
            let size_u32: u32 = a_size.max(b_size) as u32;
            if size_u32 == 1 {
                0
            } else {
                calculate_max_dna_clustering_distance(
                    size_u32,
                    sequencing_error,
                    max_clustering_distance
                )
            }
        } else {
            max_clustering_distance
        }
    };

    let a_pos1: ReferencePosition = a.get_position_1();
    let a_pos2: ReferencePosition = a.get_position_2();
    let b_pos1: ReferencePosition = b.get_position_1();
    let b_pos2: ReferencePosition = b.get_position_2();

    let distance_1: u32 = a_pos1.abs_diff(b_pos1);
    let distance_2: u32 = a_pos2.abs_diff(b_pos2);

    if distance_1 > max_distance || distance_2 > max_distance {
        return false;
    }

    if *a_variant_type == VariantType::Insertion {
        // The normalized edit distance must be within the allowed limit if they are both insertions.
        let a_sequence: String = a.get_standardized_sequence();
        let b_sequence: String = b.get_standardized_sequence();
        let max_len: f64 = a_sequence.len().max(b_sequence.len()) as f64;
        let edit_distance: f64 = edit_distance(a_sequence.as_str(), b_sequence.as_str()) as f64;
        let normalized_edit_distance: f64 = edit_distance / max_len;
        if normalized_edit_distance > max_ins_norm_edit_distance {
            return false;
        }
    }

    // A substitution is the bases it spells: another base at the same position is another
    // allele, not a second read of this one.
    if matches!(a_variant_type, VariantType::SingleNucleotideVariant | VariantType::MultiNucleotideVariant)
        && a.get_standardized_sequence() != b.get_standardized_sequence() {
        return false;
    }

    // If none of the above conditions was met, then the two variant records can be clustered.
    true
}


pub(crate) fn rebuild_orphan_pool_as_insertions(
    variant_records: &[Arc<VariantRecord>],
    indices: &[usize]
) -> Vec<Arc<VariantRecord>> {
    indices
        .iter()
        .map(|&i| {
            let vr: &VariantRecord = variant_records[i].as_ref();
            let op: &GraphOperation = vr.get_graph_operation();
            let anchor: ReferencePosition = op.get_position_1();
            let (position_1, position_2): (ReferencePosition, ReferencePosition) = match op.get_operation_type_1() {
                GraphOperationType::Upstream => (anchor.saturating_sub(1), anchor),
                _ => (anchor, anchor + 1)
            };
            Arc::new(VariantRecord::new(
                vr.get_read_id(),
                vr.get_read_position_1(),
                vr.get_read_position_2(),
                GraphOperation::new(
                    op.get_chromosome_1(),
                    position_1,
                    op.get_strand_1().clone(),
                    GraphOperationType::Downstream,
                    op.get_chromosome_1(),
                    position_2,
                    op.get_strand_1().clone(),
                    GraphOperationType::Upstream,
                    op.get_sequence().into(),
                    VariantType::Insertion
                )
            ))
        })
        .collect()
}


fn split_variant_records(
    variant_records: Vec<Arc<VariantRecord>>,
    use_thread_pool: bool,
    num_threads: usize
) -> HashMap<(ReferenceChromosomeID, ReferenceChromosomeID, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> {
    // Helper function.
    let key_of = |variant_record: &Arc<VariantRecord>| (
        variant_record.get_chromosome_1(),
        variant_record.get_chromosome_2(),
        variant_record.get_variant_type().clone(),
        variant_record.get_operation_1().clone(),
        variant_record.get_operation_2().clone()
    );

    // Step 1. Split variant records by chromosome.
    let mut variant_records_map: HashMap<(ReferenceChromosomeID, ReferenceChromosomeID, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> =
        if use_thread_pool == false {
            let thread_pool = rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .unwrap();
            thread_pool.install(|| {
                variant_records
                    .par_iter()
                    .map(|variant_record| (key_of(variant_record), Arc::clone(variant_record)))
                    .fold(
                        || HashMap::new(),
                        |mut acc, (key, variant_record)| {
                            acc.entry(key).or_insert_with(Vec::new).push(variant_record);
                            acc
                        }
                    )
                    .reduce(
                        || HashMap::new(),
                        |mut map1, map2| {
                            for (key, mut vec) in map2 {
                                map1.entry(key).or_insert_with(Vec::new).append(&mut vec);
                            }
                            map1
                        }
                    )
            })
        } else {
            let mut map: HashMap<(ReferenceChromosomeID, ReferenceChromosomeID, VariantType, GraphOperationType, GraphOperationType), Vec<Arc<VariantRecord>>> = HashMap::new();
            for variant_record in variant_records.iter() {
                map.entry(key_of(variant_record))
                    .or_insert_with(Vec::new)
                    .push(Arc::clone(variant_record));
            }
            map
        };

    // Step 2. Sort variant records by position_1.
    for records in variant_records_map.values_mut() {
        records.sort_by(|variant_record_1, variant_record_2| {
            variant_record_1.get_position_1().cmp(&variant_record_2.get_position_1())
        });
    }

    variant_records_map
}


#[cfg(test)]
#[path = "../../tests/calling/dna/variant_record_clustering.rs"]
mod tests;