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
use csv::ReaderBuilder;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use serde::de::DeserializeOwned;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


fn load_records<T: DeserializeOwned>(tsv_file: &str) -> Vec<T> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .unwrap_or_else(|error| panic!("Failed to open TSV file {}: {}", tsv_file, error));
    reader
        .deserialize()
        .map(|result| result.unwrap_or_else(|error| panic!("Failed to deserialize a row of {}: {}", tsv_file, error)))
        .collect()
}


pub fn load_rna_read_cluster_id_records(tsv_file: &str) -> Vec<RNAReadClusterIDRecord> {
    load_records(tsv_file)
}


pub fn load_rna_read_cluster_reference_transcript_records(tsv_file: &str) -> Vec<RNAReadClusterReferenceTranscriptRecord> {
    load_records(tsv_file)
}


pub fn load_rna_read_cluster_splice_junction_records(tsv_file: &str) -> Vec<RNAReadClusterSpliceJunctionRecord> {
    load_records(tsv_file)
}


pub fn load_rna_read_cluster_variant_records(tsv_file: &str) -> Vec<RNAReadClusterVariantRecord> {
    load_records(tsv_file)
}


pub fn load_rna_read_cluster_set(
    cluster_id_tsv_file: &str,
    reference_transcripts_tsv_file: &str,
    splice_junctions_tsv_file: &str,
    variants_tsv_file: &str
) -> RNAReadClusterSet {
    // Step 1. Load the four record tables.
    let cluster_id_records: Vec<RNAReadClusterIDRecord> =
        load_rna_read_cluster_id_records(cluster_id_tsv_file);
    let reference_transcript_records: Vec<RNAReadClusterReferenceTranscriptRecord> =
        load_rna_read_cluster_reference_transcript_records(reference_transcripts_tsv_file);
    let splice_junction_records: Vec<RNAReadClusterSpliceJunctionRecord> =
        load_rna_read_cluster_splice_junction_records(splice_junctions_tsv_file);
    let variant_records: Vec<RNAReadClusterVariantRecord> =
        load_rna_read_cluster_variant_records(variants_tsv_file);

    // Step 2. Rebuild the chromosome name <-> id map from every chromosome the splice-junction
    //         and variant tables name. BTreeSet keeps the assignment deterministic.
    let mut chromosome_names: BTreeSet<Box<str>> = BTreeSet::new();
    for record in splice_junction_records.iter() {
        chromosome_names.insert(record.chromosome_1.clone());
        chromosome_names.insert(record.chromosome_2.clone());
    }
    for record in variant_records.iter() {
        chromosome_names.insert(record.chromosome_1.clone());
        chromosome_names.insert(record.chromosome_2.clone());
    }
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    for (chromosome_id, chromosome_name) in chromosome_names.into_iter().enumerate() {
        chromosome_names_map.insert(chromosome_name, chromosome_id as u16);
    }
    let chromosome_id = |chromosome_name: &Box<str>| -> u16 {
        *chromosome_names_map.get_by_left(chromosome_name).unwrap()
    };

    // Step 3. Rebuild the read name <-> id map from the cluster-id table.
    let read_names: BTreeSet<Box<str>> = cluster_id_records
        .iter()
        .map(|record| record.read_name.clone())
        .collect();
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.into_iter().enumerate() {
        read_names_map.insert(read_name, read_id);
    }

    // Step 4. Group each table by cluster id. Row order within a cluster is preserved: the
    //         splice junctions form a chain in genomic order.
    let mut read_ids_by_cluster: HashMap<usize, HashSet<usize>> = HashMap::new();
    let mut shared_read_ids_by_cluster: HashMap<usize, HashSet<usize>> = HashMap::new();
    for record in cluster_id_records.iter() {
        let read_id: usize = *read_names_map.get_by_left(&record.read_name).unwrap();
        read_ids_by_cluster.entry(record.cluster_id).or_default().insert(read_id);
        if record.is_shared {
            shared_read_ids_by_cluster.entry(record.cluster_id).or_default().insert(read_id);
        }
    }
    let mut reference_gene_transcript_ids_by_cluster: HashMap<usize, HashSet<(Box<str>, Box<str>)>> = HashMap::new();
    for record in reference_transcript_records.iter() {
        reference_gene_transcript_ids_by_cluster
            .entry(record.cluster_id)
            .or_default()
            .insert((record.reference_gene_id.clone(), record.reference_transcript_id.clone()));
    }
    let mut splice_junctions_by_cluster: HashMap<usize, Vec<SpliceJunction>> = HashMap::new();
    for record in splice_junction_records.iter() {
        splice_junctions_by_cluster
            .entry(record.cluster_id)
            .or_default()
            .push(SpliceJunction::new(
                chromosome_id(&record.chromosome_1),
                chromosome_id(&record.chromosome_2),
                record.position_1,
                record.position_2,
                Strand::from_str(&record.strand_1).expect("Failed to parse strand_1"),
                Strand::from_str(&record.strand_2).expect("Failed to parse strand_2")
            ));
    }
    let mut variants_by_cluster: HashMap<usize, Vec<GraphOperation>> = HashMap::new();
    for record in variant_records.iter() {
        variants_by_cluster
            .entry(record.cluster_id)
            .or_default()
            .push(GraphOperation::new(
                chromosome_id(&record.chromosome_1),
                record.position_1,
                Strand::from_str(&record.strand_1).expect("Failed to parse strand_1"),
                GraphOperationType::from_str(&record.operation_type_1).expect("Failed to parse operation_type_1"),
                chromosome_id(&record.chromosome_2),
                record.position_2,
                Strand::from_str(&record.strand_2).expect("Failed to parse strand_2"),
                GraphOperationType::from_str(&record.operation_type_2).expect("Failed to parse operation_type_2"),
                record.sequence.clone(),
                VariantType::from_str(&record.variant_type).expect("Failed to parse variant_type")
            ));
    }

    // Step 5. Rebuild the clusters in ascending cluster-id order.
    let mut cluster_ids: Vec<usize> = read_ids_by_cluster.keys().copied().collect();
    cluster_ids.sort_unstable();
    let scores = crate::options::RNAVariantCallingOptions::DEFAULT;
    let mut variant_call_id: usize = 0;
    let mut cluster_set: RNAReadClusterSet = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    for cluster_id in cluster_ids.iter() {
        let read_ids: HashSet<usize> = read_ids_by_cluster.remove(cluster_id).unwrap();
        let witness: usize = *read_ids.iter().min().unwrap();
        let variant_calls: Vec<VariantCall> = variants_by_cluster
            .remove(cluster_id)
            .unwrap_or_default()
            .into_iter()
            .map(|graph_operation| {
                let variant_call: VariantCall = VariantCall::from_variant_records(
                    variant_call_id,
                    HashSet::from([VariantRecord::new(witness, 0, 0, graph_operation)]),
                    scores.poa_match_score,
                    scores.poa_mismatch_score,
                    scores.poa_gap_open_score,
                    scores.poa_gap_extend_score
                );
                variant_call_id += 1;
                variant_call
            })
            .collect();
        cluster_set.add_cluster(RNAReadCluster::new(
            *cluster_id,
            read_ids,
            splice_junctions_by_cluster.remove(cluster_id).unwrap_or_default(),
            variant_calls,
            HashMap::new(),
            reference_gene_transcript_ids_by_cluster.remove(cluster_id).unwrap_or_default()
        ).with_shared_read_ids(shared_read_ids_by_cluster.remove(cluster_id).unwrap_or_default()));
    }

    // Step 6. Anything still left keys on a cluster with no reads.
    assert!(
        reference_gene_transcript_ids_by_cluster.is_empty(),
        "{} references unknown cluster IDs: {:?}",
        reference_transcripts_tsv_file,
        reference_gene_transcript_ids_by_cluster.keys().collect::<Vec<_>>()
    );
    assert!(
        splice_junctions_by_cluster.is_empty(),
        "{} references unknown cluster IDs: {:?}",
        splice_junctions_tsv_file,
        splice_junctions_by_cluster.keys().collect::<Vec<_>>()
    );
    assert!(
        variants_by_cluster.is_empty(),
        "{} references unknown cluster IDs: {:?}",
        variants_tsv_file,
        variants_by_cluster.keys().collect::<Vec<_>>()
    );

    cluster_set
}


#[cfg(test)]
#[path = "../tests/io/loaders.rs"]
mod tests;
