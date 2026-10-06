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
use exacto_core::prelude::{overlaps, UnionFind};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ops::Range;

use crate::reference::splice_junction_annotation_index::SpliceJunctionAnnotationIndex;


pub(crate) fn group_transcript_models_by_shared_junctions<'a>(
    summaries: impl IntoIterator<Item = &'a RNAReadCharacterizationSummary>,
    annotation_index: &SpliceJunctionAnnotationIndex,
    max_locus_gap: u32,
    unspliced_bin_size: u32,
    min_reads: usize,
    max_slippage_repeat_length: u32,
    expected_sequencing_error: f64,
    expected_slippage_rate: f64,
    max_fpr: f64,
    num_threads: usize
) -> Vec<Vec<&'a RNAReadCharacterizationSummary>> {
    // Step 0. Annotated gene spans per chromosome, sorted, the same loci Step 9 of the clusterer
    // keys junction-less reads by, so a single-exon gene's reads load together.
    // HashMap<chromosome ID, Vec<(gene ID, start, end)>>
    let gene_spans: HashMap<u16, Vec<(Box<str>, u32, u32)>> = annotation_index
        .gene_spans
        .iter()
        .map(|(chromosome_id, spans)| {
            let mut spans: Vec<(Box<str>, u32, u32)> = spans
                .iter()
                .map(|(gene_id, &(start, end))| (gene_id.clone(), start, end))
                .collect();
            spans.sort();
            (*chromosome_id, spans)
        })
        .collect();

    // Step 1. One union-find node per distinct junction, per annotated gene, and per
    // (chromosome, bin). A spliced read unions its junctions within each locus and anchors on
    // its largest locus. A junction-less read unions the genes it overlaps and anchors on the
    // first; with no gene, it anchors on its start bin.
    let mut junction_nodes: HashMap<(u16, u32, u32), u32> = HashMap::new();
    let mut gene_nodes: HashMap<(u16, Box<str>), u32> = HashMap::new();
    let mut bin_nodes: HashMap<(u16, u32), u32> = HashMap::new();
    let mut next_node: u32 = 0;
    let mut union_find: UnionFind = UnionFind::new();
    let mut anchors: Vec<(&'a RNAReadCharacterizationSummary, u32)> = Vec::new();
    // The node of every locus of every spliced read.
    let mut read_locus_nodes: Vec<Vec<u32>> = Vec::new();
    for summary in summaries {
        let anchor: u32 = if summary.splice_junctions.is_empty() {
            let Some(chromosome) = summary.chromosome else {
                continue;
            };
            let (start, end): (u32, u32) = (
                summary.reference_start.min(summary.reference_end),
                summary.reference_start.max(summary.reference_end)
            );
            let mut anchor: Option<u32> = None;
            if let Some(spans) = gene_spans.get(&chromosome) {
                for (gene_id, gene_start, gene_end) in spans.iter() {
                    if overlaps(*gene_start as isize, *gene_end as isize, start as isize, end as isize) {
                        let node: u32 = *gene_nodes.entry((chromosome, gene_id.clone())).or_insert_with(|| {
                            next_node += 1;
                            next_node - 1
                        });
                        match anchor {
                            None => {
                                union_find.union(node, node);
                                anchor = Some(node);
                            }
                            Some(anchor) => union_find.union(anchor, node)
                        }
                    }
                }
            }
            anchor.unwrap_or_else(|| {
                let node: u32 = *bin_nodes.entry((chromosome, start / unspliced_bin_size)).or_insert_with(|| {
                    next_node += 1;
                    next_node - 1
                });
                union_find.union(node, node);
                node
            })
        } else {
            // Loci: junctions in genomic order, split at a chromosome change or a gap wider
            // than `max_locus_gap` between one intron's end and the next intron's start.
            let mut keys: Vec<(u16, u32, u32)> = summary
                .splice_junctions
                .iter()
                .map(|junction| junction.intron_span_key())
                .collect();
            keys.sort_unstable();
            let mut loci: Vec<Vec<(u16, u32, u32)>> = vec![vec![keys[0]]];
            for key in keys.into_iter().skip(1) {
                let (chromosome, _, high): (u16, u32, u32) = *loci.last().unwrap().last().unwrap();
                if key.0 == chromosome && key.1.saturating_sub(high) <= max_locus_gap {
                    loci.last_mut().unwrap().push(key);
                } else {
                    loci.push(vec![key]);
                }
            }
            // The read lives with the locus carrying most of its structure; ties go to
            // genomic order, so the choice is deterministic.
            let largest: usize = (0..loci.len())
                .max_by_key(|&i| (loci[i].len(), std::cmp::Reverse(i)))
                .unwrap();
            let mut anchor: Option<u32> = None;
            let mut locus_nodes: Vec<u32> = Vec::with_capacity(loci.len());
            for (i, locus) in loci.iter().enumerate() {
                let mut locus_anchor: Option<u32> = None;
                for key in locus.iter() {
                    let node: u32 = *junction_nodes.entry(*key).or_insert_with(|| {
                        next_node += 1;
                        next_node - 1
                    });
                    match locus_anchor {
                        None => {
                            union_find.union(node, node);
                            locus_anchor = Some(node);
                        }
                        Some(locus_anchor) => union_find.union(locus_anchor, node)
                    }
                }
                locus_nodes.push(locus_anchor.unwrap());
                if i == largest {
                    anchor = locus_anchor;
                }
            }
            read_locus_nodes.push(locus_nodes);
            anchor.unwrap()
        };
        anchors.push((summary, anchor));
    }

    // Step 1-1. Union each unannotated junction with the annotated introns it overlaps, when both
    // of its ends lie in a gene holding one of them.
    let mut annotated_introns: HashMap<u16, Vec<(u32, u32)>> = HashMap::new();
    for &(chromosome_id, start, end) in junction_nodes.keys() {
        if annotation_index.junctions_to_genes.contains_key(&(chromosome_id, start, end)) {
            annotated_introns.entry(chromosome_id).or_default().push((start, end));
        }
    }
    let mut max_intron_ends: HashMap<u16, Vec<u32>> = HashMap::new();
    for (chromosome_id, introns) in annotated_introns.iter_mut() {
        introns.sort_unstable();
        let ends: &mut Vec<u32> = max_intron_ends.entry(*chromosome_id).or_default();
        for &(_, end) in introns.iter() {
            ends.push(ends.last().map_or(end, |&max_end| max_end.max(end)));
        }
    }
    let mut novel_junctions: Vec<((u16, u32, u32), u32)> = junction_nodes
        .iter()
        .filter(|(key, _)| !annotation_index.junctions_to_genes.contains_key(*key))
        .map(|(key, node)| (*key, *node))
        .collect();
    novel_junctions.sort_unstable();
    for ((chromosome_id, start, end), node) in novel_junctions {
        let (Some(introns), Some(ends)) = (annotated_introns.get(&chromosome_id), max_intron_ends.get(&chromosome_id)) else {
            continue;
        };
        let last: usize = introns.partition_point(|&(intron_start, _)| intron_start <= end);
        for index in (0..last).rev() {
            if ends[index] < start {
                break;
            }
            let (intron_start, intron_end): (u32, u32) = introns[index];
            if intron_end < start {
                continue;
            }
            let is_within_gene: bool = annotation_index.junctions_to_genes[&(chromosome_id, intron_start, intron_end)]
                .iter()
                .any(|gene_id| {
                    annotation_index
                        .gene_spans
                        .get(&chromosome_id)
                        .and_then(|spans| spans.get(gene_id))
                        .is_some_and(|&(gene_start, gene_end)| gene_start <= start && end <= gene_end)
                });
            if is_within_gene {
                union_find.union(node, junction_nodes[&(chromosome_id, intron_start, intron_end)]);
            }
        }
    }

    // Step 2. Union the loci that enough reads bridge.
    // HashMap<node, locus>
    let mut locus_of_node: HashMap<u32, u32> = HashMap::new();
    for (locus, members) in union_find.get_clusters().into_iter().enumerate() {
        for node in members {
            locus_of_node.insert(node, locus as u32);
        }
    }
    // HashMap<locus, (node of the locus, number of reads with a junction in it)>
    let mut locus_reads: HashMap<u32, (u32, u32)> = HashMap::new();
    // BTreeMap<(locus, locus), number of reads with a junction in both>
    let mut num_bridging_reads: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for locus_nodes in read_locus_nodes.iter() {
        let mut loci: BTreeSet<u32> = BTreeSet::new();
        for node in locus_nodes.iter() {
            let locus: u32 = locus_of_node[node];
            locus_reads.entry(locus).or_insert((*node, 0));
            loci.insert(locus);
        }
        for locus in loci.iter() {
            locus_reads.get_mut(locus).unwrap().1 += 1;
        }
        for (i, locus_1) in loci.iter().enumerate() {
            for locus_2 in loci.iter().skip(i + 1) {
                *num_bridging_reads.entry((*locus_1, *locus_2)).or_insert(0) += 1;
            }
        }
    }
    if !num_bridging_reads.is_empty() {
        let depths: HashSet<u32> = locus_reads
            .values()
            .map(|(_, num_reads)| *num_reads)
            .collect();
        let min_read_support_index: RNAVariantReadSupportIndex = RNAVariantReadSupportIndex::new(
            &depths,
            max_slippage_repeat_length,
            expected_sequencing_error,
            expected_slippage_rate,
            max_fpr,
            num_threads
        );
        for ((locus_1, locus_2), num_reads) in num_bridging_reads.into_iter() {
            let (node_1, depth_1): (u32, u32) = locus_reads[&locus_1];
            let (node_2, depth_2): (u32, u32) = locus_reads[&locus_2];
            let min_read_support: u32 = min_read_support_index
                .get_min_read_support((0, depth_1.max(depth_2)))
                .max(min_reads as u32);
            if num_reads >= min_read_support {
                union_find.union(node_1, node_2);
            }
        }
    }

    // Step 3. Component per node, then every read to its anchor's component.
    let mut component_of_node: Vec<u32> = vec![0; next_node as usize];
    let mut groups: Vec<Vec<&'a RNAReadCharacterizationSummary>> = Vec::new();
    for (component, members) in union_find.get_clusters().into_iter().enumerate() {
        for node in members {
            component_of_node[node as usize] = component as u32;
        }
        groups.push(Vec::new());
    }
    for (summary, anchor) in anchors.into_iter() {
        groups[component_of_node[anchor as usize] as usize].push(summary);
    }
    groups.retain(|group| !group.is_empty());
    for group in groups.iter_mut() {
        group.sort_by_key(|summary| summary.read_id);
    }
    groups.sort_by_key(|group| group[0].read_id);
    groups
}


pub(crate) fn pack_batches(summary_groups: &[Vec<&RNAReadCharacterizationSummary>], max_batch_reads: usize) -> Vec<Range<usize>> {
    let mut batches: Vec<Range<usize>> = Vec::new();
    let (mut start, mut reads): (usize, usize) = (0, 0);
    for (i, group) in summary_groups.iter().enumerate() {
        if reads > 0 && reads + group.len() > max_batch_reads {
            batches.push(start..i);
            (start, reads) = (i, 0);
        }
        reads += group.len();
    }
    if start < summary_groups.len() {
        batches.push(start..summary_groups.len());
    }
    batches
}


#[cfg(test)]
#[path = "../tests/clustering/grouping.rs"]
mod tests;
