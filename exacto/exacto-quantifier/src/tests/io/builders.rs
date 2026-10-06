
use std::collections::{HashMap, HashSet};

use bimap::BiMap;
use exacto_caller::prelude::{GraphOperation, GraphOperationType, VariantCall, VariantRecord, VariantType};
use exacto_cluster::prelude::{RNAReadCluster, RNAReadClusterSet};
use exacto_core::prelude::Strand;

use super::*;


fn read_set(names: &[&str]) -> HashSet<Box<str>> {
    names.iter().map(|name| Box::<str>::from(*name)).collect()
}


fn pairs(items: &[(&str, &str)]) -> HashSet<(Box<str>, Box<str>)> {
    items
        .iter()
        .map(|(gene, transcript)| (Box::<str>::from(*gene), Box::<str>::from(*transcript)))
        .collect()
}


fn empty_cluster_set() -> RNAReadClusterSet {
    RNAReadClusterSet::new(BiMap::new(), BiMap::new())
}


/// The builders read gene/transcript ids and variant calls off the cluster; read ids, splice
/// junctions and genotypes are irrelevant to them, so fixtures leave those empty. Each variant
/// becomes a one-record call: with a single record the consensus is that record's own
/// operation, so the alignment scores never come into play.
fn add_cluster(
    cluster_set: &mut RNAReadClusterSet,
    id: usize,
    gene_transcript_pairs: &[(&str, &str)],
    variants: Vec<GraphOperation>
) {
    let variant_calls: Vec<VariantCall> = variants
        .into_iter()
        .enumerate()
        .map(|(variant_id, variant)| VariantCall::from_variant_records(
            variant_id,
            HashSet::from([VariantRecord::new(0, 0, 0, variant)]),
            2, -4, 4, 2
        ))
        .collect();
    cluster_set.add_cluster(RNAReadCluster::new(
        id,
        HashSet::new(),
        Vec::new(),
        variant_calls,
        HashMap::new(),
        pairs(gene_transcript_pairs)
    ));
}


/// An arbitrary, valid variant used only to make `get_variants()` non-empty.
fn single_nucleotide_variant() -> GraphOperation {
    GraphOperation::new(
        1,
        100,
        Strand::Forward,
        GraphOperationType::Mark,
        1,
        100,
        Strand::Forward,
        GraphOperationType::Mark,
        "A".into(),
        VariantType::SingleNucleotideVariant
    )
}


#[test]
fn cluster_annotation_as_str() {
    assert_eq!(ClusterAnnotation::Reference.as_str(), "reference");
    assert_eq!(ClusterAnnotation::Novel.as_str(), "novel");
    assert_eq!(ClusterAnnotation::Variant.as_str(), "variant");
}


/// Clusters 1 and 2 share a reference transcript and must roll up into one row
/// whose counts and CPM are the sums of both.
#[test]
fn build_reference_transcript_records_groups_and_sums() {
    let mut quantification_set = ClusterQuantificationSet::new();
    quantification_set.add(ClusterQuantification::new(1, 300_000.0, read_set(&["r1", "r2"])));
    quantification_set.add(ClusterQuantification::new(2, 200_000.0, read_set(&["r3"])));
    quantification_set.add(ClusterQuantification::new(3, 500_000.0, read_set(&["r4"])));

    let mut cluster_set = empty_cluster_set();
    add_cluster(&mut cluster_set, 1, &[("ENSGA", "ENSTA")], Vec::new());
    add_cluster(&mut cluster_set, 2, &[("ENSGA", "ENSTA")], Vec::new()); // same key as 1
    add_cluster(&mut cluster_set, 3, &[("ENSGB", "ENSTB")], Vec::new());

    let records: Vec<QuantificationReferenceTranscriptRecord> =
        build_quantification_reference_transcript_records(&quantification_set, &cluster_set)
            .collect();

    assert_eq!(records.len(), 2);

    let group_a = records.iter().find(|r| &*r.reference_gene_id == "ENSGA").unwrap();
    assert_eq!(group_a.num_clusters, 2);
    assert_eq!(group_a.num_reads, 3); // 2 + 1
    assert!((group_a.cpm - 500_000.0).abs() < 1e-9); // 300k + 200k

    let group_b = records.iter().find(|r| &*r.reference_gene_id == "ENSGB").unwrap();
    assert_eq!(group_b.num_clusters, 1);
    assert_eq!(group_b.num_reads, 1);
    assert!((group_b.cpm - 500_000.0).abs() < 1e-9);
}


/// Two clusters carry the same gene/transcript pair set in opposite order; sorting
/// the pairs must collapse them onto one reference-transcript row (canonical key).
#[test]
fn build_reference_transcript_records_uses_canonical_key() {
    let mut quantification_set = ClusterQuantificationSet::new();
    quantification_set.add(ClusterQuantification::new(1, 100_000.0, read_set(&["r1"])));
    quantification_set.add(ClusterQuantification::new(2, 100_000.0, read_set(&["r2"])));

    let mut cluster_set = empty_cluster_set();
    add_cluster(&mut cluster_set, 1, &[("G1", "T1"), ("G2", "T2")], Vec::new());
    add_cluster(&mut cluster_set, 2, &[("G2", "T2"), ("G1", "T1")], Vec::new());

    let records: Vec<QuantificationReferenceTranscriptRecord> =
        build_quantification_reference_transcript_records(&quantification_set, &cluster_set)
            .collect();

    assert_eq!(records.len(), 1, "reversed pair orders must map to one key");
    assert_eq!(&*records[0].reference_gene_id, "G1;G2");
    assert_eq!(&*records[0].reference_transcript_id, "T1;T2");
    assert_eq!(records[0].num_clusters, 2);
    assert!((records[0].cpm - 200_000.0).abs() < 1e-9);
}


/// Exercises every gene/transcript/allele derivation branch:
///   1: reference gene + reference transcript, reference allele
///   2: reference gene + novel transcript (empty transcript id), reference allele
///   3: novel gene (no reference overlap), novel transcript, reference allele
///   4: reference gene + reference transcript, variant allele (carries a variant)
#[test]
fn build_cluster_annotation_records_derives_all_states() {
    let mut quantification_set = ClusterQuantificationSet::new();
    for id in 1..=4 {
        quantification_set.add(ClusterQuantification::new(id, 250_000.0, read_set(&["r"])));
    }

    let mut cluster_set = empty_cluster_set();
    add_cluster(&mut cluster_set, 1, &[("ENSG1", "ENST1")], Vec::new());
    add_cluster(&mut cluster_set, 2, &[("ENSG2", "")], Vec::new());
    add_cluster(&mut cluster_set, 3, &[], Vec::new());
    add_cluster(&mut cluster_set, 4, &[("ENSG4", "ENST4")], vec![single_nucleotide_variant()]);

    let records: Vec<QuantificationClusterRecord> =
        build_quantification_cluster_records(&quantification_set, &cluster_set)
            .collect();
    let by_id: HashMap<usize, &QuantificationClusterRecord> =
        records.iter().map(|record| (record.cluster_id, record)).collect();

    let c1 = by_id[&1];
    assert_eq!(&*c1.gene_annotation, "reference");
    assert_eq!(&*c1.transcript_annotation, "reference");
    assert_eq!(&*c1.allele, "reference");
    assert_eq!(&*c1.reference_gene_id, "ENSG1");
    assert_eq!(&*c1.reference_transcript_id, "ENST1");

    let c2 = by_id[&2];
    assert_eq!(&*c2.gene_annotation, "reference");
    assert_eq!(&*c2.transcript_annotation, "novel");
    assert_eq!(&*c2.allele, "reference");
    assert_eq!(&*c2.reference_gene_id, "ENSG2");
    assert_eq!(&*c2.reference_transcript_id, "");

    let c3 = by_id[&3];
    assert_eq!(&*c3.gene_annotation, "novel");
    assert_eq!(&*c3.transcript_annotation, "novel");
    assert_eq!(&*c3.allele, "reference");
    assert_eq!(&*c3.reference_gene_id, "");
    assert_eq!(&*c3.reference_transcript_id, "");

    let c4 = by_id[&4];
    assert_eq!(&*c4.gene_annotation, "reference");
    assert_eq!(&*c4.transcript_annotation, "reference");
    assert_eq!(&*c4.allele, "variant");
}
