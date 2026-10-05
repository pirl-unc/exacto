use bimap::BiMap;
use exacto_core::prelude::{FastaMap, Strand};
use std::fs;
use std::path::Path;

use crate::prelude::{GraphOperationType, VariantRecord};

use super::*;


#[test]
fn is_repeat_variant_deletion_returns_matches() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_837,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_841,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);

    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_840,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);

    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_836,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);

    let go: GraphOperation = GraphOperation::new(
        0,
        4_330_358,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_330_361,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, false);
}

#[test]
fn is_repeat_variant_insertion_returns_matches() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let go: GraphOperation = GraphOperation::new(
        0,
        3_829_837,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        3_829_838,
        Strand::Forward,
        GraphOperationType::Upstream,
        "AAAA".into(),
        VariantType::Insertion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        104,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 9);

    let go: GraphOperation = GraphOperation::new(
        0,
        4_329_639,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_329_640,
        Strand::Forward,
        GraphOperationType::Upstream,
        "AT".into(),
        VariantType::Insertion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 8);

    let go: GraphOperation = GraphOperation::new(
        0,
        4_329_640,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        4_329_641,
        Strand::Forward,
        GraphOperationType::Upstream,
        "TA".into(),
        VariantType::Insertion
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 6);
}

#[test]
fn is_repeat_variant_snv_returns_matches() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let go: GraphOperation = GraphOperation::new(
        0,
        7_674_222,
        Strand::Reverse,
        GraphOperationType::Downstream,
        0,
        7_674_225,
        Strand::Reverse,
        GraphOperationType::Upstream,
        "GG".into(),
        VariantType::MultiNucleotideVariant
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, false);
    assert_eq!(size, 0);

    let go: GraphOperation = GraphOperation::new(
        0,
        7_670_078,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        7_670_081,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CC".into(),
        VariantType::MultiNucleotideVariant
    );

    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, true);
    assert_eq!(size, 5);
}

#[test]
fn is_repeat_variant_mnv_returns_matches() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    let go: GraphOperation = GraphOperation::new(
        0,
        7_674_002,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        7_674_005,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CC".into(),
        VariantType::MultiNucleotideVariant
    );
    let vr: VariantRecord = VariantRecord::new(
        1,
        101,
        102,
        go
    );

    let (is_repeat, size) = is_repeat_variant(
        &vr.get_graph_operation(),
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, false);
    assert_eq!(size, 0);
}

/// A deletion of one base in sequence that holds no repeat is not a repeat variant.
///
///   chr17:7,674,077-7,674,085   G T A G T A G T A
///   Deleted                             T           (7,674,081)
///
/// No two adjacent bases are equal and no two adjacent dinucleotides are equal.
#[test]
fn is_repeat_variant_returns_false_for_deletion_of_one_base_outside_repeats() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    assert_eq!(fasta_map.get_sequence("chr17", 7_674_077, 7_674_085).to_uppercase(), "GTAGTAGTA");

    let go: GraphOperation = GraphOperation::new(
        0,
        7_674_080,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        7_674_082,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );

    let (is_repeat, size) = is_repeat_variant(
        &go,
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, false);
    assert_eq!(size, 0);
}

/// An insertion of one base in sequence that holds no repeat is not a repeat variant.
///
///   chr17:7,674,077-7,674,085   G T A G T A G T A
///   Inserted                             C          (between 7,674,081 and 7,674,082)
#[test]
fn is_repeat_variant_returns_false_for_insertion_of_one_base_outside_repeats() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();

    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"),0u16)].into_iter().collect();
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    assert_eq!(fasta_map.get_sequence("chr17", 7_674_077, 7_674_085).to_uppercase(), "GTAGTAGTA");

    let go: GraphOperation = GraphOperation::new(
        0,
        7_674_081,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        7_674_082,
        Strand::Forward,
        GraphOperationType::Upstream,
        "C".into(),
        VariantType::Insertion
    );

    let (is_repeat, size) = is_repeat_variant(
        &go,
        &chromosome_names_map,
        &fasta_map,
        4,
        4
    );

    assert_eq!(is_repeat, false);
    assert_eq!(size, 0);
}
