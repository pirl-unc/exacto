use exacto_annotator::prelude::*;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_variant_integration_1() {
    let tsv_path_1 = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants_annotated.tsv");
    let tsv_full_path_1 = fs::canonicalize(tsv_path_1).unwrap();
    let tsv_file_1: &str = tsv_full_path_1.to_str().unwrap();
    let tsv_path_2 = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv");
    let tsv_full_path_2 = fs::canonicalize(tsv_path_2).unwrap();
    let tsv_file_2: &str = tsv_full_path_2.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();
    let dna_variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(tsv_file_1);
    let rna_variant_call_set: RNAVariantCallSet = RNAVariantCallSet::read_tsv_file(tsv_file_2);
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let integrated_variant_set: IntegratedVariantSet = integrate_dna_rna_variants(
        &dna_variant_call_annotation_set,
        &rna_variant_call_set,
        &gene_annotator,
        2,
        1000,
        100_000,
        1
    );

    assert_eq!(integrated_variant_set.get_size(), 6);

    let mut rna_variant_call_ids: HashSet<usize> = HashSet::new();
    for integrated_variant in integrated_variant_set.integrated_variants.iter() {
        rna_variant_call_ids.insert(integrated_variant.rna_variant_call_id);
    }

    assert_eq!(rna_variant_call_ids.contains(&1), true);
    assert_eq!(rna_variant_call_ids.contains(&3), true);
    assert_eq!(rna_variant_call_ids.contains(&5), true);
    assert_eq!(rna_variant_call_ids.contains(&10), true);
    assert_eq!(rna_variant_call_ids.contains(&12), true);
    assert_eq!(rna_variant_call_ids.contains(&17), true);
}

