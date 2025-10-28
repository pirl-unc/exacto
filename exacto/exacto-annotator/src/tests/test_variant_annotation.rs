use exacto_core::prelude::*;
use polars::prelude::*;
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_variant_annotation_1() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
}

#[test]
fn test_variant_annotation_2() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-002-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
}

#[test]
fn test_variant_annotation_3() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-003-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.exon_ids.contains_key("ENST00000504937.5"), true);
}

#[test]
fn test_variant_annotation_4() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-004-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141510.18"), true);
}

#[test]
fn test_variant_annotation_5() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-005-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141510.18"), true);
}

#[test]
fn test_variant_annotation_6() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-006-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000231824.4"), true);
}

#[test]
fn test_variant_annotation_7() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-007-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000231824.4"), true);
}

#[test]
fn test_variant_annotation_8() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-008-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141499.18"), true);
}

#[test]
fn test_variant_annotation_9() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-009-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141499.18"), true);
}

#[test]
fn test_variant_annotation_10() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-010-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141499.18"), true);
}

#[test]
fn test_variant_annotation_11() {
    let tsv_path = Path::new("src/tests/data/tsv/dna-011-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let variant_call_annotation_set: VariantCallAnnotationSet = annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.gene_ids.contains("ENSG00000141499.18"), true);
}
