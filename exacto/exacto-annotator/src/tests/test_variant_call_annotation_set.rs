use exacto_core::prelude::*;
use polars::prelude::*;
use std::fs;
use std::path::Path;
use tempfile::NamedTempFile;

use crate::prelude::*;


#[test]
fn test_variant_call_annotation_set_1() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(tsv_file);

    assert_eq!(variant_call_annotation_set.annotations.keys().len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
}

#[test]
fn test_variant_call_annotation_set_2() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let variant_call_annotation_set_1: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(tsv_file);
    let variant_call_annotation_set_2: VariantCallAnnotationSet = variant_call_annotation_set_1.clone();

    assert_eq!(variant_call_annotation_set_1, variant_call_annotation_set_2);
}

#[test]
fn test_variant_call_annotation_set_3() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let variant_call_annotation_set_1: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(tsv_file);
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    variant_call_annotation_set_1.to_tsv_file(file.path().to_str().unwrap());

    let variant_call_annotation_set_2: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set_1, variant_call_annotation_set_2);
}

#[test]
fn test_variant_call_annotation_set_4() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(tsv_file);

    assert_eq!(variant_call_annotation_set.get(1).variant_call_id, 1);
    assert_eq!(variant_call_annotation_set.get_by_range("chr17", 7_000_000, 8_000_000).len(), 1);
    assert_eq!(variant_call_annotation_set.get_by_reference_transcript("ENST00000269305.9").is_some(), true);
}
