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


use std::any::Any;
use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use exacto_core::log_info;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


pub fn build_genome_variation_graph(
    fasta_file: &str,
    df_variants: &DataFrame,
    graph_type: VarGraphTypes,
    num_threads: usize
) -> Vec<VarGraph> {
    // Step 1. Get DataFrame columns
    let col_variant_id = df_variants.column("variant_id").unwrap().i64().unwrap();
    let col_chromosome_1 = df_variants.column("chromosome_1").unwrap().str().unwrap();
    let col_position_1 = df_variants.column("position_1").unwrap().i64().unwrap();
    let col_operation_1 = df_variants.column("operation_1").unwrap().str().unwrap();
    let col_strand_1 = df_variants.column("strand_1").unwrap().str().unwrap();
    let col_chromosome_2 = df_variants.column("chromosome_2").unwrap().str().unwrap();
    let col_position_2 = df_variants.column("position_2").unwrap().i64().unwrap();
    let col_operation_2 = df_variants.column("operation_2").unwrap().str().unwrap();
    let col_strand_2 = df_variants.column("strand_2").unwrap().str().unwrap();
    let col_sequence = df_variants.column("sequence").unwrap().str().unwrap();
    let col_variant_type = df_variants.column("variant_type").unwrap().str().unwrap();

    let num_cycles_series = df_variants
        .column("num_cycles")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();

    let col_num_cycles = num_cycles_series.str().unwrap();

    // Step 2. Create GraphOperation objects
    log_info!("Started loading variants.");
    let mut graph_operations_map: HashMap<usize, GraphOperationView> = HashMap::new();
    for i in 0..df_variants.height() {
        let variant_id: usize = col_variant_id.get(i).unwrap() as usize;
        let chromosome_1: &str = col_chromosome_1.get(i).unwrap();
        let position_1: usize = col_position_1.get(i).unwrap() as usize;
        let operation_type_1: GraphOperationType = GraphOperationType::from_str(col_operation_1.get(i).unwrap()).unwrap();
        let strand_1: Strand = Strand::from_str(col_strand_1.get(i).unwrap()).unwrap();
        let chromosome_2: &str = col_chromosome_2.get(i).unwrap();
        let position_2: usize = col_position_2.get(i).unwrap() as usize;
        let operation_type_2: GraphOperationType = GraphOperationType::from_str(col_operation_2.get(i).unwrap()).unwrap();
        let strand_2: Strand = Strand::from_str(col_strand_2.get(i).unwrap()).unwrap();
        let sequence: &str = col_sequence.get(i).unwrap_or("");
        let num_cycles: Option<u16> = col_num_cycles
            .get(i)
            .and_then(|s| {
                let s = s.trim();
                if s.is_empty() {
                    return None;
                }

                // Try integer first
                if let Ok(v) = s.parse::<u16>() {
                    return Some(v);
                }

                // Try float next (e.g. "2.0")
                if let Ok(v) = s.parse::<f64>() {
                    if v.is_finite() && v.fract() == 0.0 && v >= 0.0 {
                        return Some(v as u16);
                    }
                }

                None
            });

        let gov: GraphOperationView = GraphOperationView::new(
            &*chromosome_1,
            position_1 as u32,
            &operation_type_1,
            &strand_1,
            &*chromosome_2,
            position_2 as u32,
            &operation_type_2,
            &strand_2,
            &*sequence,
            num_cycles
        );
        graph_operations_map.insert(variant_id, gov);
    }
    log_info!("Finished loading variants.");

    // Step 3. Build a map between chromosome names and IDs
    let mut chromosome_lengths_map: HashMap<Box<str>, u32> = HashMap::new();
    let mut chromosome_ids_map: BiMap<Box<str>, usize> = BiMap::new();
    for (i, (chromosome, length)) in get_fasta_sequence_ids(fasta_file).iter().enumerate() {
        chromosome_lengths_map.insert(chromosome.clone(), *length);
        chromosome_ids_map.insert(chromosome.clone(), i);
    }

    // Step 4. Cluster reference chromosomes
    log_info!("Started clustering reference chromosomes.");
    let mut uf: UnionFind = UnionFind::new();
    for (chromosome, chromosome_id) in chromosome_ids_map.iter() {
        uf.union(*chromosome_id as u32, *chromosome_id as u32);
    }
    for (variant_id, graph_op_view) in graph_operations_map.iter() {
        let chromosome_1_id: usize = *chromosome_ids_map.get_by_left(graph_op_view.get_chromosome_1()).unwrap();
        let chromosome_2_id: usize = *chromosome_ids_map.get_by_left(graph_op_view.get_chromosome_2()).unwrap();
        uf.union(chromosome_1_id as u32, chromosome_2_id as u32);
    }
    log_info!("Finished clustering reference chromosomes.");

    // Step 5. Assign a variation graph ID for each cluster of reference chromosomes
    let mut vargraph_chromosomes_map: HashMap<usize, Vec<Box<str>>> = HashMap::new();
    let mut vargraph_ids_map: HashMap<Box<str>, usize> = HashMap::new();
    let mut vargraph_id: usize = 0;
    for cluster in uf.get_clusters().iter() {
        for chromosome_id in cluster.iter() {
            let chromosome: Box<str> = chromosome_ids_map.get_by_right(&(*chromosome_id as usize)).unwrap().clone();
            vargraph_chromosomes_map
                .entry(vargraph_id)
                .or_insert(Vec::new())
                .push(chromosome.clone());
            vargraph_ids_map.insert(chromosome.clone(), vargraph_id);
        }
        vargraph_id += 1;
    }
    for (chromosome, chromosome_id) in chromosome_ids_map.iter() {
        if vargraph_ids_map.contains_key(chromosome) == false {
            vargraph_chromosomes_map
                .entry(vargraph_id)
                .or_insert(Vec::new())
                .push(chromosome.clone());
            vargraph_id += 1;
        }
    }

    // Step 6. Assign each GraphOperationView to the corresponding variation graph ID
    log_info!("Started assigning each GraphOperationView to the corresponding variation graph ID.");
    let mut vargraph_operations_map: HashMap<usize, (HashSet<Box<str>>, HashMap<usize, GraphOperationView>)> = HashMap::new();
    for (vargraph_id, chromosomes) in vargraph_chromosomes_map.iter() {
        vargraph_operations_map.insert(*vargraph_id, (chromosomes.iter().cloned().collect(), HashMap::new()));
    }
    for (variant_id, graph_op_view) in graph_operations_map.iter() {
        let chromosome_1: Box<str> = graph_op_view.get_chromosome_1().into();
        let chromosome_2: Box<str> = graph_op_view.get_chromosome_2().into();
        let vargraph_id_1: usize = *vargraph_ids_map.get(&chromosome_1).unwrap();
        let vargraph_id_2: usize = *vargraph_ids_map.get(&chromosome_2).unwrap();
        assert_eq!(vargraph_id_1, vargraph_id_2);
        vargraph_operations_map.get_mut(&vargraph_id_1).unwrap().1.insert(*variant_id, graph_op_view.clone());
    }
    log_info!("Started assigning each GraphOperationView to the corresponding variation graph ID.");

    // Step 7. Add variants to the variation graphs
    log_info!("Started adding variants to the variation graphs.");
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut vargraphs: Vec<VarGraph> = thread_pool.install(|| {
        vargraph_operations_map
            .par_iter()
            .map(|(vargraph_id, (reference_chromosomes, graph_op_views_map))| {
                let reference_chromosomes_string: String = reference_chromosomes
                    .iter()
                    .map(|s| s.as_ref())
                    .collect::<Vec<&str>>()
                    .join(",");

                // Collect variant breakpoints (positions where we want 1-bp reference nodes)
                let mut breakpoint_positions: HashMap<Box<str>, BTreeSet<u32>> = HashMap::new();
                for (_variant_id, gov) in graph_op_views_map.iter() {
                    let chromosome_1: &str = gov.get_chromosome_1();
                    let position_1: u32 = gov.get_position_1();
                    breakpoint_positions
                        .entry(chromosome_1.into())
                        .or_insert_with(BTreeSet::new)
                        .insert(position_1);

                    let chromosome_2: &str = gov.get_chromosome_2();
                    let position_2: u32 = gov.get_position_2();
                    breakpoint_positions
                        .entry(chromosome_2.into())
                        .or_insert_with(BTreeSet::new)
                        .insert(position_2);
                }

                // Build pre-split VarGraphReferenceNode lists for each chromosome
                let mut reference_nodes: Vec<Vec<VarGraphReferenceNode>> = Vec::new();
                for chromosome in reference_chromosomes.iter() {
                    let chr: &str = chromosome.as_ref();
                    let chr_len: u32 = *chromosome_lengths_map
                        .get(chromosome)
                        .unwrap_or_else(|| panic!("Missing chromosome length for {chr}"));

                    // Full sequence for this chromosome
                    let seq: Box<str> = get_fasta_sequence(chr, 1, chr_len, fasta_file);

                    // Breakpoints for this chromosome, sorted and deduped
                    let mut breakpoints: Vec<u32> = breakpoint_positions
                        .get(chr)
                        .map(|set| set.iter().copied().collect())
                        .unwrap_or_else(Vec::new);
                    breakpoints.sort_unstable();
                    breakpoints.dedup();

                    let mut chromosome_reference_nodes: Vec<VarGraphReferenceNode> = Vec::new();
                    let mut curr_pos: u32 = 1;
                    let mut seq_offset: usize = 0;
                    for bp in breakpoints {
                        // Region [curr_pos .. bp-1], if any
                        if curr_pos < bp {
                            let size: u32 = bp - curr_pos;
                            let end: u32 = bp - 1;
                            let slice_end: usize = seq_offset + size as usize;
                            let subseq: &str = &seq[seq_offset..slice_end];
                            chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                                chr,
                                curr_pos,
                                end,
                                subseq,
                            ));
                            seq_offset = slice_end;
                            curr_pos = bp;
                        }

                        // Single base [bp .. bp]
                        if curr_pos == bp {
                            let slice_end: usize = seq_offset + 1;
                            let subseq: &str = &seq[seq_offset..slice_end];
                            chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                                chr,
                                bp,
                                bp,
                                subseq,
                            ));
                            seq_offset = slice_end;
                            curr_pos = bp + 1;
                        }
                    }

                    // Tail region [curr_pos .. chr_len], if any
                    if curr_pos <= chr_len {
                        let size: u32 = chr_len - curr_pos + 1;
                        let slice_end: usize = seq_offset + size as usize;
                        let subseq: &str = &seq[seq_offset..slice_end];
                        chromosome_reference_nodes.push(VarGraphReferenceNode::new(
                            chr,
                            curr_pos,
                            chr_len,
                            subseq
                        ));
                    }
                    reference_nodes.push(chromosome_reference_nodes);
                }

                let mut vargraph: VarGraph = VarGraph::from_reference_nodes(reference_nodes);

                // Add variants
                let mut count: usize = 0;
                let total: usize = graph_op_views_map.len();
                for (variant_id, graph_op_view) in graph_op_views_map.iter() {
                    count += 1;
                    if count % 10_000 == 0 {
                        log_info!("VarGraph for {reference_chromosomes_string} processed {count}/{total}.");
                    }
                    let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(
                        *variant_id,
                        graph_op_view.clone()
                    );
                    vargraph.add_variant_node(variant_node);
                }

                if graph_type == VarGraphTypes::Individual {
                    // Only stitch adjacent variants if we're building an individual (personalized) graph
                    vargraph.stitch_adjacent_variants();
                }

                vargraph
            })
            .collect()
    });
    log_info!("Finished adding variants to the variation graphs.");

    vargraphs
}

/// Builds a single transcriptome variation graph from a pre-partitioned group DataFrame.
pub fn build_single_transcriptome_variation_graph(
    fasta_file: &str,
    df_group: &DataFrame,
    chromosome_lengths_map: &HashMap<Box<str>, u32>,
    graph_type: VarGraphTypes,
) -> (usize, Box<str>, VarGraph) {
    let col_transcript_model_id = df_group
        .column("transcript_model_id")
        .unwrap()
        .i64()
        .unwrap();
    let col_index = df_group.column("index").unwrap().i64().unwrap();
    let col_sequence = df_group.column("sequence").unwrap().str().unwrap();
    let col_chromosome_1 = df_group.column("chromosome_1").unwrap().str().unwrap();
    let col_position_1 = df_group.column("position_1").unwrap().i64().unwrap();
    let col_operation_1 = df_group.column("operation_1").unwrap().str().unwrap();
    let col_strand_1 = df_group.column("strand_1").unwrap().str().unwrap();
    let col_chromosome_2 = df_group.column("chromosome_2").unwrap().str().unwrap();
    let col_position_2 = df_group.column("position_2").unwrap().i64().unwrap();
    let col_operation_2 = df_group.column("operation_2").unwrap().str().unwrap();
    let col_strand_2 = df_group.column("strand_2").unwrap().str().unwrap();

    let num_cycles_series = df_group
        .column("num_cycles")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();

    let col_num_cycles = num_cycles_series.str().unwrap();

    let transcript_model_id: usize = col_transcript_model_id.get(0).unwrap() as usize;

    let reference_transcript_ids: Box<str> = df_group
        .column("reference_transcript_ids")
        .unwrap()
        .cast(&DataType::String)
        .unwrap()
        .str()
        .unwrap()
        .get(0)
        .unwrap_or("")
        .into();

    // Get GraphOperationViews from the DataFrame
    let mut graph_op_views: Vec<(usize, GraphOperationView)> = Vec::new();
    for i in 0..df_group.height() {
        let index: usize = col_index.get(i).unwrap() as usize;
        let chromosome_1: &str = col_chromosome_1.get(i).unwrap();
        let position_1: u32 = col_position_1.get(i).unwrap() as u32;
        let operation_type_1: GraphOperationType = GraphOperationType::from_str(col_operation_1.get(i).unwrap()).unwrap();
        let strand_1: Strand = Strand::from_str(col_strand_1.get(i).unwrap()).unwrap();
        let chromosome_2: &str = col_chromosome_2.get(i).unwrap();
        let position_2: u32 = col_position_2.get(i).unwrap() as u32;
        let operation_type_2: GraphOperationType = GraphOperationType::from_str(col_operation_2.get(i).unwrap()).unwrap();
        let strand_2: Strand = Strand::from_str(col_strand_2.get(i).unwrap()).unwrap();
        let num_cycles: Option<u16> = col_num_cycles
            .get(i)
            .and_then(|v| v.parse::<u16>().ok());
        let mut sequence: Box<str> = col_sequence.get(i).unwrap_or("").into();

        if operation_type_1 == GraphOperationType::Include &&
            operation_type_2 == GraphOperationType::Include &&
            sequence == "".into() {
            sequence = get_fasta_sequence(
                chromosome_1,
                position_1,
                position_2,
                fasta_file
            );
        }

        let graph_op_view: GraphOperationView = GraphOperationView::new(
            chromosome_1,
            position_1,
            &operation_type_1,
            &strand_1,
            chromosome_2,
            position_2,
            &operation_type_2,
            &strand_2,
            &*sequence,
            num_cycles
        );

        assert!(i == index);

        graph_op_views.push((index, graph_op_view));
    }

    graph_op_views.sort_by_key(|(idx, _)| *idx);

    // Build per-chromosome (start, end) spans from graph operations
    let mut chromosome_spans: HashMap<Box<str>, (u32, u32)> = HashMap::new();
    for (_, gov) in graph_op_views.iter() {
        for (chr, pos) in [
            (gov.get_chromosome_1(), gov.get_position_1()),
            (gov.get_chromosome_2(), gov.get_position_2()),
        ] {
            let entry = chromosome_spans
                .entry(chr.into())
                .or_insert((pos, pos));
            entry.0 = entry.0.min(pos);
            entry.1 = entry.1.max(pos);
        }
    }

    // Load only the needed region per chromosome
    let mut reference_nodes: Vec<Vec<VarGraphReferenceNode>> = Vec::new();
    for (chromosome, (min_pos, max_pos)) in chromosome_spans.iter() {
        let chromosome_length: u32 = *chromosome_lengths_map.get(chromosome.as_ref()).unwrap();
        let start: u32 = (*min_pos).max(1);
        let end: u32   = (*max_pos).min(chromosome_length);
        let sequence: Box<str> = get_fasta_sequence(chromosome, start, end, fasta_file);
        let reference_node = VarGraphReferenceNode::new(chromosome, start, end, &*sequence);
        reference_nodes.push(vec![reference_node]);
    }

    // Create a VarGraph from the reference nodes
    let mut vargraph: VarGraph = VarGraph::from_reference_nodes(reference_nodes);

    // Add variants
    for (variant_id, graph_op_view) in graph_op_views.iter() {
        let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(
            *variant_id,
            graph_op_view.clone()
        );
        vargraph.add_variant_node(variant_node);
    }

    // Stitch adjacent variants
    if graph_type == VarGraphTypes::Individual {
        // Only stitch adjacent variants if we're building an individual (personalized) graph
        vargraph.stitch_adjacent_variants();
    }

    (transcript_model_id, reference_transcript_ids, vargraph)
}

pub fn build_transcriptome_variation_graph(
    fasta_file: &str,
    df_transcript_structures: &DataFrame,
    graph_type: VarGraphTypes,
    num_threads: usize
) -> Vec<(usize, Box<str>, VarGraph)> {
    let mut chromosome_lengths_map: HashMap<Box<str>, u32> = HashMap::new();
    for (chromosome, length) in get_fasta_sequence_ids(fasta_file).iter() {
        chromosome_lengths_map.insert(chromosome.clone(), *length);
    }
    let groups: Vec<DataFrame> = df_transcript_structures
        .partition_by(["transcript_model_id", "reference_transcript_ids"], true)
        .unwrap();
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = thread_pool.install(|| {
        groups
            .par_iter()
            .map(|df_group| {
                build_single_transcriptome_variation_graph(
                    fasta_file,
                    df_group,
                    &chromosome_lengths_map,
                    graph_type.clone()
                )
            })
            .collect()
    });

    vargraphs
}
