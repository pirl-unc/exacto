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
use exacto_core::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use tempfile::TempPath;

use crate::clustering::breakpoint_clip::resolve_breakpoint_clips;
use crate::clustering::splice_junction_cluster::SpliceJunctionCluster;
use crate::clustering::splice_junction_clustering::cluster_transcript_models_by_splice_junctions;
use crate::options::ClusterRNAReadsOptions;
use crate::phasing::cells::store_cells;
use crate::phasing::genotype_phasing::phase_genotypes;
use crate::read_cluster::rna_read_cluster_set::RNAReadClusterSet;
use crate::read_filtering::read_cleanup::clean_transcript_model;
use crate::reference::splice_junction_annotation_index::SpliceJunctionAnnotationIndex;
use crate::variant_filtering::cluster_depths::ClusterDepths;
use crate::variant_filtering::known_variants::KnownVariants;
use crate::variant_filtering::variant_judge::{min_reads_at_site, VariantJudge};


pub(crate) struct BatchContext<'a, A: GeneAnnotator + Sync> {
    pub(crate) gene_annotator: &'a A,
    pub(crate) chromosome_names_map: &'a BiMap<Box<str>, u16>,
    pub(crate) fasta_map: &'a FastaMap,
    pub(crate) annotation_index: &'a SpliceJunctionAnnotationIndex,
    
    /// HashMap<(chromosome ID, intron boundary), Vec<intron>> of the introns the reads splice.
    pub(crate) intron_boundary_index: &'a HashMap<(u16, u32), Vec<(u32, u32)>>,
    
    pub(crate) known_variants: &'a KnownVariants<'a>,

    /// The read support of a variant at every depth a junction cluster can hold.
    pub(crate) read_support_index: &'a RNAVariantReadSupportIndex,
    
    /// The transcript models `characterize_rna_reads` spilled, read back per summary group.
    pub(crate) models_file: &'a TempPath,
    pub(crate) options: &'a ClusterRNAReadsOptions,
    pub(crate) thread_pool: &'a ThreadPool
}


pub(crate) fn cluster_batch<A: GeneAnnotator + Sync>(
    context: &BatchContext<A>,
    summary_groups: &[Vec<&RNAReadCharacterizationSummary>],
    next_cluster_id: &mut usize,
    next_variant_call_id: &mut usize,
    rna_read_cluster_set: &mut RNAReadClusterSet
) {
    let BatchContext { gene_annotator, chromosome_names_map, fasta_map, known_variants, read_support_index, options, thread_pool, .. } = *context;

    // Load and cluster each group of the batch in parallel.
    let loaded: Vec<(Vec<TranscriptModel>, Vec<SpliceJunctionCluster>)> = thread_pool.install(|| {
        summary_groups
            .par_iter()
            .map(|summary_group| cluster_group(context, summary_group))
            .collect()
    });

    // Flatten the batch: cluster ids offset per group in group order, so every later
    // sort on cluster id reproduces today's group-by-group numbering exactly.
    let mut cluster_read_ids: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut splice_junction_clusters_map: HashMap<usize, Arc<SpliceJunctionCluster>> = HashMap::new();
    let mut transcript_models_map: HashMap<usize, Arc<TranscriptModel>> = HashMap::new();
    let mut offset: usize = 0;
    for (transcript_models, splice_junction_clusters) in loaded {
        let next_offset: usize = offset + splice_junction_clusters.len();
        for mut cluster in splice_junction_clusters {
            cluster.id += offset;
            cluster_read_ids.insert(cluster.id, cluster.read_ids.iter().copied().collect());
            splice_junction_clusters_map.insert(cluster.id, Arc::new(cluster));
        }
        transcript_models_map.extend(transcript_models.into_iter().map(|m| (m.get_read_id(), Arc::new(m))));
        offset = next_offset;
    }

    // Identify RNA variants within each junction cluster, and read each call's depth off the
    // cluster's own reads.
    // HashMap<cluster ID, Vec<(VariantCall, (depth at side 1, depth at side 2))>>
    log_info!("Identifying RNA variants within each junction cluster.");
    let variant_calls_by_cluster: HashMap<usize, Vec<(VariantCall, (u32, u32))>> =
        thread_pool.install(|| {
            cluster_read_ids
                .par_iter()
                .map(|(cluster_id, read_ids)| {
                    let splice_junction_cluster: &SpliceJunctionCluster = &splice_junction_clusters_map[cluster_id];
                    (*cluster_id, call_cluster_variants(splice_junction_cluster, read_ids, &transcript_models_map, options))
                })
                .collect()
        });

    // Filter variant calls: the allow list, DNA variant support, total depth, read support,
    // and template switch signatures (see `VariantJudge`).
    log_info!("Filtering variant calls.");
    let judge: VariantJudge<A> = VariantJudge::new(
        known_variants,
        read_support_index,
        &transcript_models_map,
        gene_annotator,
        chromosome_names_map,
        fasta_map,
        &options.filtering,
        options.calling.max_ins_norm_edit_distance,
        &options.analyte_type
    );

    // Judge every cluster's variant calls in parallel.
    // Vec<(cluster ID, Vec<(VariantCall, repeat length)>, Vec<FailedVariantCall>)>
    let mut judged: Vec<(usize, Vec<(VariantCall, Option<u32>)>, Vec<FailedVariantCall>)> = thread_pool.install(|| {
        variant_calls_by_cluster
            .into_par_iter()
            .map(|(cluster_id, cluster_variant_calls)| {
                let mut kept: Vec<(VariantCall, Option<u32>)> = Vec::new();
                let mut failed: Vec<FailedVariantCall> = Vec::new();
                for (variant_call, depths) in cluster_variant_calls {
                    match judge.judge(&variant_call, depths) {
                        Ok(repeat_length) => kept.push((variant_call, repeat_length)),
                        Err(failure) => failed.push(FailedVariantCall::new(variant_call, failure))
                    }
                }
                (cluster_id, kept, failed)
            })
            .collect()
    });

    judged.sort_by_key(|(cluster_id, _, _)| *cluster_id);
    let mut variant_calls_by_cluster: HashMap<usize, Vec<VariantCall>> = HashMap::with_capacity(judged.len());
    let mut failed_by_cluster: HashMap<usize, Vec<FailedVariantCall>> = HashMap::with_capacity(judged.len());
    // HashMap<variant call ID, repeat length its read support was judged at>
    let mut repeat_lengths: HashMap<usize, Option<u32>> = HashMap::new();
    for (cluster_id, kept, failed) in judged {
        let mut variant_calls: Vec<VariantCall> = Vec::with_capacity(kept.len());
        for (mut variant_call, repeat_length) in kept {
            // Assign global variant call IDs.
            variant_call.set_id(*next_variant_call_id);
            repeat_lengths.insert(*next_variant_call_id, repeat_length);
            *next_variant_call_id += 1;
            variant_calls.push(variant_call);
        }
        variant_calls_by_cluster.insert(cluster_id, variant_calls);
        failed_by_cluster.insert(cluster_id, failed);
    }

    // Genotype the RNA variants.
    // Vec<(cluster ID, HashMap<variant call ID, Vec<(read ID, Allele)>>)>
    log_info!("Genotyping the variant calls.");
    let variant_genotypes: Vec<(usize, HashMap<usize, Vec<(usize, Allele)>>)> = thread_pool.install(|| {
        variant_calls_by_cluster
            .par_iter()
            .map(|(&cluster_id, variant_calls)| {
                let cluster: &SpliceJunctionCluster = &splice_junction_clusters_map[&cluster_id];

                // HashMap<variant call ID, Vec<(read ID, Allele)>>
                let genotypes: HashMap<usize, Vec<(usize, Allele)>> = genotype_rna_variants(
                    &cluster.read_ids,
                    variant_calls,
                    &transcript_models_map
                );

                (cluster_id, genotypes)
            })
            .collect()
    });

    // Phase the RNA variants.
    // Vec<(cluster ID, HashMap<cell ID, HashSet<read ID>>, HashMap<cell ID, HashSet<shared read ID>>)>
    log_info!("Phasing the variant calls.");

    // A cell is held to the read support of a variant at a depth of its cluster's reads,
    // and never to less than a junction cluster is.
    let mut cluster_cell_read_ids: Vec<(usize, HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>)> = thread_pool.install(|| {
        variant_genotypes
            .par_iter()
            .map(|(cluster_id, genotypes)| {
                let read_ids: &HashSet<usize> = &cluster_read_ids[cluster_id];
                let min_read_support: usize = floor_at(read_support_index, read_ids.len() as u32);
                let (cell_read_ids, shared_read_ids): (HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>) = phase_genotypes(
                    read_ids,
                    genotypes,
                    min_read_support.max(options.junction.min_reads),
                    options.phasing.mec_max_k,
                    options.phasing.mec_num_restart,
                    options.phasing.mec_max_iter,
                    options.seed as u64
                );
                (*cluster_id, cell_read_ids, shared_read_ids)
            })
            .collect()
    });
    
    // Store the results: one RNAReadCluster per phased cell.
    log_info!("Storing the clustering results.");
    let mut genotypes_by_cluster: HashMap<usize, HashMap<usize, Vec<(usize, Allele)>>> = variant_genotypes.into_iter().collect();
    cluster_cell_read_ids.sort_by_key(|(cluster_id, _, _)| *cluster_id);
    for (cluster_id, cell_read_ids, shared_read_ids) in cluster_cell_read_ids.into_iter() {
        store_cells(
            &splice_junction_clusters_map[&cluster_id],
            &variant_calls_by_cluster[&cluster_id],
            genotypes_by_cluster.remove(&cluster_id).unwrap_or_default(),
            failed_by_cluster.remove(&cluster_id).unwrap_or_default(),
            (cell_read_ids, shared_read_ids),
            &repeat_lengths,
            &transcript_models_map,
            read_support_index,
            next_cluster_id,
            rna_read_cluster_set
        );
    }
}


fn cluster_group<A: GeneAnnotator + Sync>(
    context: &BatchContext<A>,
    summary_group: &[&RNAReadCharacterizationSummary]
) -> (Vec<TranscriptModel>, Vec<SpliceJunctionCluster>) {
    let BatchContext { chromosome_names_map, fasta_map, annotation_index, intron_boundary_index, models_file, options, .. } = *context;
    let floor_at = |depth: u32| -> usize { floor_at(context.read_support_index, depth) };
    let mut transcript_models: Vec<TranscriptModel> = load_transcript_models(models_file, summary_group.iter().copied());
    for transcript_model in transcript_models.iter_mut() {
        clean_transcript_model(transcript_model, intron_boundary_index, chromosome_names_map, fasta_map, &options.calling);
    }
    
    let mut splice_junction_clusters: Vec<SpliceJunctionCluster> = cluster_transcript_models_by_splice_junctions(
        &transcript_models,
        annotation_index,
        chromosome_names_map,
        fasta_map,
        options.junction.min_reads,
        options.error_model.sequencing_error,
        options.filtering.max_fpr,
        options.junction.max_unspliced_locus_gap,
        None
    );

    if options.calling.soft_clip_removal {
        resolve_breakpoint_clips(&mut transcript_models, &mut splice_junction_clusters,
            chromosome_names_map, fasta_map, &options.calling, &options.filtering, floor_at);
        splice_junction_clusters.retain(|cluster| cluster.read_ids.len() >= options.junction.min_reads);
        for (id, cluster) in splice_junction_clusters.iter_mut().enumerate() {
            cluster.id = id;
        }
    }
    
    (transcript_models, splice_junction_clusters)
}


fn call_cluster_variants(
    splice_junction_cluster: &SpliceJunctionCluster,
    read_ids: &HashSet<usize>,
    transcript_models_map: &HashMap<usize, Arc<TranscriptModel>>,
    options: &ClusterRNAReadsOptions
) -> Vec<(VariantCall, (u32, u32))> {
    let junctions: &Vec<SpliceJunction> = splice_junction_cluster.splice_junctions.as_ref();
    let transcript_models: Vec<Arc<TranscriptModel>> = read_ids
        .iter()
        .map(|read_id| Arc::clone(transcript_models_map.get(read_id).unwrap()))
        .collect();

    // Cluster variant records among the transcript models.
    let mut cluster_variant_calls: Vec<VariantCall> = cluster_rna_variant_records(
        junctions,
        &transcript_models,
        options.calling.max_clustering_distance,
        options.calling.min_size_proportion,
        options.calling.max_ins_norm_edit_distance,
        options.error_model.sequencing_error,
        options.calling.poa_match_score,
        options.calling.poa_mismatch_score,
        options.calling.poa_gap_open_score,
        options.calling.poa_gap_extend_score
    );
    
    cluster_variant_calls.retain(|variant_call| matches!(
        variant_call.get_consensus_graph_operation().get_variant_type(),
        VariantType::SingleNucleotideVariant
            | VariantType::MultiNucleotideVariant
            | VariantType::Insertion
            | VariantType::Deletion
            | VariantType::Breakpoint
            | VariantType::CircularRNA
            | VariantType::FusionGene
            | VariantType::Translocation
    ));

    let cluster_depths: ClusterDepths = ClusterDepths::new(&transcript_models);
    let cluster_variant_calls: Vec<(VariantCall, (u32, u32))> = cluster_variant_calls
        .into_iter()
        .map(|mut variant_call| {
            let op: &GraphOperation = variant_call.get_consensus_graph_operation();
            let carriers: Vec<&TranscriptModel> = variant_call
                .get_read_ids()
                .iter()
                .map(|read_id| &**transcript_models_map.get(read_id).unwrap())
                .collect();
            let depths: (u32, u32) = cluster_depths.at_site(op, &carriers);
            variant_call.set_total_depth(depths.0.max(depths.1) as i32);
            (variant_call, depths)
        })
        .collect();

    cluster_variant_calls
}


fn floor_at(read_support_index: &RNAVariantReadSupportIndex, depth: u32) -> usize {
    min_reads_at_site(read_support_index, 0, (depth, depth)) as usize
}
