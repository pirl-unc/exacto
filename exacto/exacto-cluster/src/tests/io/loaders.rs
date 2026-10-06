use std::collections::{HashMap, HashSet};
use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use tempfile::{tempdir, TempDir};

use super::*;
use crate::io::builders::*;


/// Writing a set's four tables and loading them back gives the same clusters: the same reads,
/// reference ids, junctions and variant spellings per cluster. Chromosomes are compared by name,
/// because the two sets number them independently.
#[test]
fn test_cluster_set_round_trips_through_its_tables() {
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0);
    chromosome_names_map.insert("chr2".into(), 1);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in ["read_a", "read_b", "read_c"].iter().enumerate() {
        read_names_map.insert((*read_name).into(), read_id);
    }
    let substitution: GraphOperation = GraphOperation::new(
        0, 149, Strand::Forward, GraphOperationType::Downstream,
        0, 151, Strand::Forward, GraphOperationType::Upstream,
        "T".into(), VariantType::SingleNucleotideVariant
    );
    let breakpoint: GraphOperation = GraphOperation::new(
        0, 300, Strand::Forward, GraphOperationType::Downstream,
        1, 900, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let call = |id: usize, read_id: usize, graph_operation: GraphOperation| -> VariantCall {
        VariantCall::from_variant_records(
            id, HashSet::from([VariantRecord::new(read_id, 0, 0, graph_operation)]), 0, 4, 6, 2
        )
    };
    let mut cluster_set: RNAReadClusterSet =
        RNAReadClusterSet::new(read_names_map.clone(), chromosome_names_map.clone());
    cluster_set.add_cluster(RNAReadCluster::new(
        0,
        HashSet::from([0, 1]),
        vec![SpliceJunction::new(0, 0, 100, 200, Strand::Forward, Strand::Forward)],
        vec![call(0, 0, substitution)],
        HashMap::new(),
        HashSet::from([("GENE1".into(), "TX1".into())])
    ));
    cluster_set.add_cluster(RNAReadCluster::new(
        1,
        HashSet::from([2]),
        vec![SpliceJunction::new(0, 1, 300, 900, Strand::Forward, Strand::Reverse)],
        vec![call(1, 2, breakpoint)],
        HashMap::new(),
        HashSet::new()
    ));

    let directory: TempDir = tempdir().unwrap();
    let path = |name: &str| -> String { directory.path().join(name).to_str().unwrap().to_string() };
    write_tsv_file(build_rna_read_cluster_id_records(&cluster_set), directory.path().join("clusters.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_reference_gene_transcript_records(&cluster_set), directory.path().join("reference_transcripts.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_splice_junction_records(&cluster_set), directory.path().join("splice_junctions.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_variant_records(&cluster_set), directory.path().join("variants.tsv").as_path()).unwrap();

    let loaded: RNAReadClusterSet = load_rna_read_cluster_set(
        &path("clusters.tsv"),
        &path("reference_transcripts.tsv"),
        &path("splice_junctions.tsv"),
        &path("variants.tsv")
    );

    assert_eq!(loaded.get_clusters().len(), 2);
    let name = |set: &RNAReadClusterSet, chromosome_id: u16| -> Box<str> {
        set.get_chromosome_names_map().get_by_right(&chromosome_id).unwrap().clone()
    };
    for cluster in cluster_set.get_clusters() {
        let reloaded: &RNAReadCluster = loaded.get_cluster(cluster.get_id());
        assert_eq!(
            reloaded.get_read_names(loaded.get_read_names_map()),
            cluster.get_read_names(cluster_set.get_read_names_map())
        );
        assert_eq!(reloaded.get_reference_gene_transcript_ids(), cluster.get_reference_gene_transcript_ids());

        let junctions = |set: &RNAReadClusterSet, cluster: &RNAReadCluster| -> Vec<(Box<str>, Box<str>, u32, u32, String, String)> {
            cluster.get_splice_junctions().iter().map(|junction| (
                name(set, junction.chromosome_1),
                name(set, junction.chromosome_2),
                junction.position_1,
                junction.position_2,
                junction.strand_1.as_str().to_string(),
                junction.strand_2.as_str().to_string()
            )).collect()
        };
        assert_eq!(junctions(&loaded, reloaded), junctions(&cluster_set, cluster));

        let variants = |set: &RNAReadClusterSet, cluster: &RNAReadCluster| -> Vec<(Box<str>, u32, String, Box<str>, u32, String, String, String)> {
            cluster.get_variant_calls().iter().map(|variant_call| {
                let graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
                (
                    name(set, graph_operation.get_chromosome_1()),
                    graph_operation.get_position_1(),
                    graph_operation.get_operation_type_1().as_str().to_string(),
                    name(set, graph_operation.get_chromosome_2()),
                    graph_operation.get_position_2(),
                    graph_operation.get_operation_type_2().as_str().to_string(),
                    graph_operation.get_sequence().to_string(),
                    graph_operation.get_variant_type().as_str().to_string()
                )
            }).collect()
        };
        assert_eq!(variants(&loaded, reloaded), variants(&cluster_set, cluster));
    }
}


/// A read shared between two clusters is a row of each, marked `is_shared` in the cluster it
/// fits as well, and loads back as a shared read of that cluster only. `read_b` was assigned to
/// cluster 0 and is shared with cluster 1. A table written before the column existed loads with
/// no shared reads.
#[test]
fn test_cluster_set_round_trips_shared_reads_through_its_tables() {
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in ["read_a", "read_b", "read_c"].iter().enumerate() {
        read_names_map.insert((*read_name).into(), read_id);
    }
    let junction: SpliceJunction = SpliceJunction::new(0, 0, 100, 200, Strand::Forward, Strand::Forward);
    let mut cluster_set: RNAReadClusterSet = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(0, HashSet::from([0, 1]), vec![junction.clone()], Vec::new(), HashMap::new(), HashSet::new()));
    cluster_set.add_cluster(
        RNAReadCluster::new(1, HashSet::from([1, 2]), vec![junction], Vec::new(), HashMap::new(), HashSet::new())
            .with_shared_read_ids(HashSet::from([1]))
    );

    let rows: Vec<(usize, Box<str>, bool)> = build_rna_read_cluster_id_records(&cluster_set)
        .map(|record| (record.cluster_id, record.read_name, record.is_shared))
        .collect();
    assert_eq!(
        rows,
        vec![
            (0, "read_a".into(), false),
            (0, "read_b".into(), false),
            (1, "read_b".into(), true),
            (1, "read_c".into(), false)
        ]
    );

    let directory: TempDir = tempdir().unwrap();
    let path = |name: &str| -> String { directory.path().join(name).to_str().unwrap().to_string() };
    write_tsv_file(build_rna_read_cluster_id_records(&cluster_set), directory.path().join("clusters.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_reference_gene_transcript_records(&cluster_set), directory.path().join("reference_transcripts.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_splice_junction_records(&cluster_set), directory.path().join("splice_junctions.tsv").as_path()).unwrap();
    write_tsv_file(build_rna_read_cluster_variant_records(&cluster_set), directory.path().join("variants.tsv").as_path()).unwrap();
    std::fs::write(path("clusters_without_column.tsv"), "cluster_id\tread_name\n0\tread_a\n0\tread_b\n1\tread_c\n").unwrap();

    let loaded: RNAReadClusterSet = load_rna_read_cluster_set(
        &path("clusters.tsv"),
        &path("reference_transcripts.tsv"),
        &path("splice_junctions.tsv"),
        &path("variants.tsv")
    );
    let shared_names = |set: &RNAReadClusterSet, cluster_id: usize| -> HashSet<Box<str>> {
        set.get_cluster(cluster_id)
            .get_shared_read_ids()
            .iter()
            .map(|read_id| set.get_read_names_map().get_by_right(read_id).unwrap().clone())
            .collect()
    };
    assert_eq!(shared_names(&loaded, 0), HashSet::new());
    assert_eq!(shared_names(&loaded, 1), HashSet::from(["read_b".into()]));
    assert_eq!(loaded.get_cluster(1).get_read_ids().len(), 2);
    assert_eq!(
        loaded.get_read_clusters_map(),
        HashMap::from([("read_a".into(), 0), ("read_b".into(), 0), ("read_c".into(), 1)])
    );

    let old: RNAReadClusterSet = load_rna_read_cluster_set(
        &path("clusters_without_column.tsv"),
        &path("reference_transcripts.tsv"),
        &path("splice_junctions.tsv"),
        &path("variants.tsv")
    );
    assert!(old.get_clusters().iter().all(|cluster| cluster.get_shared_read_ids().is_empty()));
}
