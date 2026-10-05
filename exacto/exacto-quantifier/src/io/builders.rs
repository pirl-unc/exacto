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


use exacto_core::prelude::*;
use exacto_cluster::prelude::*;
use std::collections::HashMap;

use crate::prelude::*;


pub fn build_quantification_cluster_records<'a>(
    quantification_set: &'a ClusterQuantificationSet,
    rna_read_cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = QuantificationClusterRecord> + 'a {
    assert!(
        !quantification_set.quantifications.is_empty(),
        "quantification_set.quantifications is empty."
    );
    assert!(
        !rna_read_cluster_set.get_clusters().is_empty(),
        "rna_read_cluster_set.clusters is empty."
    );

    quantification_set.quantifications.iter().map(move |quantification| {
        let cluster: &RNAReadCluster = rna_read_cluster_set.get_cluster(quantification.get_cluster_id());

        // Sort so the joined id lists are canonical (mirrors the reference-transcript roll-up).
        let mut pairs: Vec<&(Box<str>, Box<str>)> = cluster
            .get_reference_gene_transcript_ids()
            .iter()
            .collect();
        pairs.sort_unstable();

        // Novel when the cluster overlaps no reference gene
        let gene_annotation: ClusterAnnotation = if pairs.iter().any(|(reference_gene_id, _)| !reference_gene_id.is_empty()) {
            ClusterAnnotation::Reference
        } else {
            ClusterAnnotation::Novel
        };

        // Novel when no pair carries a reference transcript id
        let transcript_annotation: ClusterAnnotation = if pairs.iter().any(|(_, reference_transcript_id)| !reference_transcript_id.is_empty()) {
            ClusterAnnotation::Reference
        } else {
            ClusterAnnotation::Novel
        };

        // Variant when the cluster carries at least one variant, else reference
        let allele: ClusterAnnotation = if cluster.get_variant_calls().is_empty() {
            ClusterAnnotation::Reference
        } else {
            ClusterAnnotation::Variant
        };

        let reference_gene_id: Box<str> = pairs
            .iter()
            .map(|(reference_gene_id, _)| &**reference_gene_id)
            .collect::<Vec<&str>>()
            .join(LIST_SEPARATOR)
            .into();

        let reference_transcript_id: Box<str> = pairs
            .iter()
            .map(|(_, reference_transcript_id)| &**reference_transcript_id)
            .collect::<Vec<&str>>()
            .join(LIST_SEPARATOR)
            .into();

        QuantificationClusterRecord {
            cluster_id: quantification.get_cluster_id(),
            gene_annotation: gene_annotation.as_str().into(),
            reference_gene_id,
            transcript_annotation: transcript_annotation.as_str().into(),
            reference_transcript_id,
            allele: allele.as_str().into(),
            num_reads: quantification.get_read_count(),
            cpm: quantification.get_cpm()
        }
    })
}


pub fn build_quantification_reference_transcript_records<'a>(
    quantification_set: &'a ClusterQuantificationSet,
    rna_read_cluster_set: &'a RNAReadClusterSet
) -> impl Iterator<Item = QuantificationReferenceTranscriptRecord> + 'a {
    assert!(
        !quantification_set.quantifications.is_empty(),
        "quantification_set.quantifications is empty."
    );
    assert!(
        !rna_read_cluster_set.get_clusters().is_empty(),
        "rna_read_cluster_set.clusters is empty."
    );

    // HashMap<(reference gene IDs, reference transcript IDs), (num clusters, num reads, sum of CPM)>
    let mut reference_transcript_map: HashMap<(Box<str>, Box<str>), (usize, usize, f64)> = HashMap::new();

    for quantification in quantification_set.quantifications.iter() {
        let cluster = rna_read_cluster_set.get_cluster(quantification.get_cluster_id());

        // Sort the (gene, transcript) pairs so the joined key is canonical: two
        // clusters carrying the same pair set in a different order must collapse
        // onto one reference-transcript row, not split their CPM across two keys.
        let mut pairs: Vec<&(Box<str>, Box<str>)> =
            cluster.get_reference_gene_transcript_ids().iter().collect();
        pairs.sort_unstable();

        let reference_gene_ids: Box<str> = pairs
            .iter()
            .map(|(reference_gene_id, _)| &**reference_gene_id)
            .collect::<Vec<&str>>()
            .join(LIST_SEPARATOR)
            .into();
        let reference_transcript_ids: Box<str> = pairs
            .iter()
            .map(|(_, reference_transcript_id)| &**reference_transcript_id)
            .collect::<Vec<&str>>()
            .join(LIST_SEPARATOR)
            .into();

        let entry = reference_transcript_map
            .entry((reference_gene_ids, reference_transcript_ids))
            .or_insert((0, 0, 0.0));

        entry.0 += 1;
        entry.1 += quantification.get_read_count();
        entry.2 += quantification.get_cpm();
    }

    // Sorted before emission: `reference_transcript_map` is a `HashMap`, so raw
    // iteration order would reorder this file's rows on every run even though the
    // aggregated values are identical.
    let mut rows: Vec<((Box<str>, Box<str>), (usize, usize, f64))> =
        reference_transcript_map.into_iter().collect();
    rows.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    rows.into_iter().map(
        |((reference_gene_id, reference_transcript_id), (num_clusters, num_reads, cpm))| {
            QuantificationReferenceTranscriptRecord {
                reference_gene_id,
                reference_transcript_id,
                num_clusters,
                num_reads,
                cpm
            }
        }
    )
}


#[cfg(test)]
#[path = "../tests/io/builders.rs"]
mod tests;