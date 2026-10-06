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

use crate::prelude::*;


fn sorted_clusters(cluster_set: &RNAReadClusterSet) -> Vec<&RNAReadCluster> {
    let mut clusters: Vec<&RNAReadCluster> = cluster_set.get_clusters();
    clusters.sort_by_key(|cluster| cluster.get_id());
    clusters
}


fn variant_record(
    cluster_set: &RNAReadClusterSet,
    cluster_id: usize,
    graph_operation: &GraphOperation
) -> RNAReadClusterVariantRecord {
    let chromosome_name = |chromosome_id: u16| -> Box<str> {
        cluster_set.get_chromosome_names_map().get_by_right(&chromosome_id).unwrap().clone()
    };
    RNAReadClusterVariantRecord {
        cluster_id,
        chromosome_1: chromosome_name(graph_operation.get_chromosome_1()),
        position_1: graph_operation.get_position_1(),
        strand_1: graph_operation.get_strand_1().as_str().into(),
        operation_type_1: graph_operation.get_operation_type_1().as_str().into(),
        chromosome_2: chromosome_name(graph_operation.get_chromosome_2()),
        position_2: graph_operation.get_position_2(),
        strand_2: graph_operation.get_strand_2().as_str().into(),
        operation_type_2: graph_operation.get_operation_type_2().as_str().into(),
        sequence: graph_operation.get_sequence().into(),
        variant_type: graph_operation.get_variant_type().as_str().into()
    }
}


pub fn build_rna_read_cluster_id_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterIDRecord> + 'a {
    sorted_clusters(cluster_set).into_iter().flat_map(move |cluster| {
        // Sorted, not just collected: the names come out of a set, so emitting them in
        // iteration order would make identical inputs produce different files.
        let mut read_names: Vec<(Box<str>, bool)> = cluster
            .get_read_ids()
            .iter()
            .map(|read_id| {
                let read_name: Box<str> = cluster_set.get_read_names_map().get_by_right(read_id).unwrap().clone();
                (read_name, cluster.get_shared_read_ids().contains(read_id))
            })
            .collect();
        read_names.sort_unstable();
        read_names.into_iter().map(move |(read_name, is_shared)| RNAReadClusterIDRecord {
            cluster_id: cluster.get_id(),
            read_name,
            is_shared
        })
    })
}


pub fn build_rna_read_cluster_summary_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterSummaryRecord> + 'a {
    sorted_clusters(cluster_set).into_iter().map(|cluster| RNAReadClusterSummaryRecord {
        cluster_id: cluster.get_id(),
        num_reads: cluster.get_read_ids().len()
    })
}


pub fn build_rna_read_cluster_reference_gene_transcript_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterReferenceTranscriptRecord> + 'a {
    sorted_clusters(cluster_set).into_iter().flat_map(|cluster| {
        let mut pairs: Vec<&(Box<str>, Box<str>)> = cluster.get_reference_gene_transcript_ids().iter().collect();
        pairs.sort_unstable();
        pairs.into_iter().map(move |(reference_gene_id, reference_transcript_id)| {
            RNAReadClusterReferenceTranscriptRecord {
                cluster_id: cluster.get_id(),
                reference_gene_id: reference_gene_id.clone(),
                reference_transcript_id: reference_transcript_id.clone()
            }
        })
    })
}


pub fn build_rna_read_cluster_splice_junction_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterSpliceJunctionRecord> + 'a {
    sorted_clusters(cluster_set).into_iter().flat_map(move |cluster| {
        cluster.get_splice_junctions().iter().map(move |splice_junction| {
            let chromosome_name = |chromosome_id: u16| -> Box<str> {
                cluster_set.get_chromosome_names_map().get_by_right(&chromosome_id).unwrap().clone()
            };
            RNAReadClusterSpliceJunctionRecord {
                cluster_id: cluster.get_id(),
                chromosome_1: chromosome_name(splice_junction.chromosome_1),
                chromosome_2: chromosome_name(splice_junction.chromosome_2),
                position_1: splice_junction.position_1,
                position_2: splice_junction.position_2,
                strand_1: splice_junction.strand_1.as_str().into(),
                strand_2: splice_junction.strand_2.as_str().into()
            }
        })
    })
}


pub fn build_rna_read_cluster_variant_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterVariantRecord> + 'a {
    sorted_clusters(cluster_set).into_iter().flat_map(move |cluster| {
        cluster.get_variant_calls().iter().map(move |variant_call| {
            variant_record(cluster_set, cluster.get_id(), variant_call.get_consensus_graph_operation())
        })
    })
}


pub fn build_rna_read_cluster_failed_variant_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterFailedVariantRecord> + 'a {
    cluster_set.get_failed_variant_calls().iter().map(move |failed| {
        let num_reads: u32 = failed.num_reads;
        let total_depth: Option<u32> = u32::try_from(failed.total_depth).ok();
        let (alt_count, coverage, slippage_repeat_length): (u32, Option<u32>, Option<u32>) = match &failed.failure {
            VariantCallFailure::LowTotalDepth { total_depth, .. } => (num_reads, Some(*total_depth), None),
            VariantCallFailure::TooFewReads { num_reads, repeat_len, .. } => (*num_reads, total_depth, *repeat_len),
            _ => (num_reads, total_depth, None)
        };
        let columns: RNAReadClusterVariantRecord = variant_record(cluster_set, failed.cluster_id, &failed.graph_operation);
        RNAReadClusterFailedVariantRecord {
            cluster_id: Some(failed.cluster_id),
            chromosome_1: columns.chromosome_1,
            position_1: columns.position_1,
            strand_1: columns.strand_1,
            operation_type_1: columns.operation_type_1,
            chromosome_2: columns.chromosome_2,
            position_2: columns.position_2,
            strand_2: columns.strand_2,
            operation_type_2: columns.operation_type_2,
            sequence: columns.sequence,
            variant_type: columns.variant_type,
            failure_reason: failed.failure.as_str().into(),
            alt_count,
            reference_count: coverage.map(|coverage| coverage.saturating_sub(alt_count)).unwrap_or(0),
            coverage,
            p_value: None,
            p_value_cutoff: None,
            expected_error_rate: None,
            slippage_repeat_length,
            junction_homology: None
        }
    })
}


pub fn build_rna_read_cluster_template_switch_records<'a>(
    cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = RNAReadClusterTemplateSwitchRecord> + 'a {
    cluster_set.get_failed_variant_calls().iter().filter_map(move |failed| {
        let VariantCallFailure::TemplateSwitch(reason) = &failed.failure else {
            return None;
        };
        let columns: RNAReadClusterVariantRecord = variant_record(cluster_set, failed.cluster_id, &failed.graph_operation);
        Some(RNAReadClusterTemplateSwitchRecord {
            cluster_id: failed.cluster_id,
            chromosome_1: columns.chromosome_1,
            position_1: columns.position_1,
            strand_1: columns.strand_1,
            operation_type_1: columns.operation_type_1,
            chromosome_2: columns.chromosome_2,
            position_2: columns.position_2,
            strand_2: columns.strand_2,
            operation_type_2: columns.operation_type_2,
            sequence: columns.sequence,
            variant_type: columns.variant_type,
            homology_left: None,
            homology_right: None,
            canonical_splice: None,
            breakpoint_dispersion: None,
            foldback: matches!(reason, TemplateSwitchReason::Foldback),
            partner_count: None,
            num_reads: failed.num_reads,
            flagged: true,
            suppressed: true
        })
    })
}


#[cfg(test)]
#[path = "../tests/io/builders.rs"]
mod tests;