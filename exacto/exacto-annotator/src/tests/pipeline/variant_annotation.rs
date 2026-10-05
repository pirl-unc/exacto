use std::collections::HashSet;
use std::fs;
use std::path::Path;
use tempfile::NamedTempFile;

use super::*;


#[test]
fn annotate_variant_calls_returns_matches_for_dna_001() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-001-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_002() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-002-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_003() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-003-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_004() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-004-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_005() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-005-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003725258.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003723991.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_006() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-006-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 2);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141499.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000698743.1").unwrap().as_ref(), "ENSE00003974626.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intergenic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141499.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000698743.1").unwrap().as_ref(), "ENSE00003974626.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intergenic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_gene_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_007() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-007-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000108381.11"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141255.13"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000179314.16"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000576947.1").unwrap().as_ref(), "ENSE00002647845.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_008() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-008-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_009() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-009-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003786593.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_010() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-010-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_011() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-011-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003786593.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_012() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-012-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000420246.6").unwrap().as_ref(), "ENSE00003735852.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_013() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-013-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 2);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000420246.6").unwrap().as_ref(), "ENSE00003735852.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_014() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-014-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00002037735.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00002037735.1");
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_015() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-015-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 2);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000108381.11"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141255.13"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000179314.16"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000576947.1").unwrap().as_ref(), "ENSE00002647845.1");
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000179314.16"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000072818.12"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&2usize).unwrap().position_2_annotation.reference_exon_ids.is_empty(), true);
}

#[test]
fn annotate_variant_calls_returns_matches_for_dna_016() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-016-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.genic_region, GenicRegion::Intronic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000108381.11"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_gene_ids.contains("ENSG00000141255.13"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_1_annotation.reference_exon_ids.is_empty(), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_gene_ids.contains("ENSG00000179314.16"), true);
    assert_eq!(variant_call_annotation_set.annotations.get(&1usize).unwrap().position_2_annotation.reference_exon_ids.get("ENST00000576947.1").unwrap().as_ref(), "ENSE00002647845.1");
}

/// The table is built by the variant caller's own writer, so this test fails when the caller
/// and the annotator stop agreeing on a column name.
#[test]
fn annotate_variant_calls_reads_table_written_by_variant_caller() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let mut df_variant_records = dna_variant_records_to_dataframe(vec![
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 7,
            chromosome_1: "chr17".into(),
            position_1: 7674224,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674226,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read_1;read_2".into(),
            num_consensus_read_names: 2,
            read_names: "read_1;read_2".into(),
            num_read_names: 2
        }
    ]);

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    let mut tsv = fs::File::create(file.path()).unwrap();
    CsvWriter::new(&mut tsv)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df_variant_records)
        .unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
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

    assert_eq!(variant_call_annotation_set.annotations.len(), 1);
    assert_eq!(variant_call_annotation_set.get(7).id, 7);
    assert_eq!(variant_call_annotation_set.get(7).position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.get(7).position_1_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
    assert_eq!(variant_call_annotation_set.to_dataframe().get_column_names()[0].as_str(), "variant_id");
}

#[test]
fn append_variant_call_annotations_keeps_columns_and_rows_in_order() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    // Variant IDs are out of order. The deletion has an empty sequence.
    let mut df_variant_records = dna_variant_records_to_dataframe(vec![
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 3,
            chromosome_1: "chr17".into(),
            position_1: 7674224,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674226,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read_1;read_2".into(),
            num_consensus_read_names: 2,
            read_names: "read_1;read_2".into(),
            num_read_names: 2
        },
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 1,
            chromosome_1: "chr17".into(),
            position_1: 7674200,
            strand_1: "-".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674231,
            strand_2: "-".into(),
            operation_2: "U".into(),
            sequence: "".into(),
            variant_size: Some(30),
            variant_type: "DEL".into(),
            consensus_read_names: "read_3".into(),
            num_consensus_read_names: 1,
            read_names: "read_3;read_4".into(),
            num_read_names: 2
        },
        DNAVariantRecord {
            origin: "germline".into(),
            variant_id: 2,
            chromosome_1: "chr17".into(),
            position_1: 7676155,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr18".into(),
            position_2: 5170100,
            strand_2: "-".into(),
            operation_2: "D".into(),
            sequence: "".into(),
            variant_size: None,
            variant_type: "TRA".into(),
            consensus_read_names: "read_5".into(),
            num_consensus_read_names: 1,
            read_names: "read_5".into(),
            num_read_names: 1
        }
    ]);

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    let mut tsv = fs::File::create(file.path()).unwrap();
    CsvWriter::new(&mut tsv)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df_variant_records)
        .unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    // Every input column, in its order, then the annotation columns in the order of the positions
    let mut expected_column_names: Vec<String> = df_variant_calls
        .get_column_names()
        .iter()
        .map(|name| name.to_string())
        .collect();
    expected_column_names.push("position_1_genic_region".to_string());
    expected_column_names.push("position_1_annotation".to_string());
    expected_column_names.push("position_1_plus_1_genic_region".to_string());
    expected_column_names.push("position_1_plus_1_annotation".to_string());
    expected_column_names.push("position_2_minus_1_genic_region".to_string());
    expected_column_names.push("position_2_minus_1_annotation".to_string());
    expected_column_names.push("position_2_genic_region".to_string());
    expected_column_names.push("position_2_annotation".to_string());
    let column_names: Vec<String> = df_annotated
        .get_column_names()
        .iter()
        .map(|name| name.to_string())
        .collect();
    assert_eq!(column_names, expected_column_names);

    // Every input cell is unchanged
    assert_eq!(df_annotated.height(), 3);
    assert_eq!(df_annotated.select(df_variant_calls.get_column_names_owned()).unwrap().equals_missing(&df_variant_calls), true);

    // Every row carries its own annotation
    let col_position_1_genic_region = df_annotated.column("position_1_genic_region").unwrap().str().unwrap();
    let col_position_2_genic_region = df_annotated.column("position_2_genic_region").unwrap().str().unwrap();
    let col_position_2_annotation = df_annotated.column("position_2_annotation").unwrap().str().unwrap();
    assert_eq!(col_position_1_genic_region.get(0).unwrap(), "exonic");
    assert_eq!(col_position_2_genic_region.get(0).unwrap(), "exonic");
    assert_eq!(col_position_1_genic_region.get(1).unwrap(), "exonic");
    assert_eq!(col_position_2_genic_region.get(1).unwrap(), "exonic");
    assert_eq!(col_position_1_genic_region.get(2).unwrap(), "exonic");
    assert_eq!(col_position_2_genic_region.get(2).unwrap(), "intronic");
    assert_eq!(col_position_2_annotation.get(0).unwrap().starts_with("ENSG00000141510.18|"), true);
    assert_eq!(col_position_2_annotation.get(2).unwrap().starts_with("ENSG00000231824.4;ENSG00000266153.1|"), true);
}

#[test]
fn append_variant_call_annotations_output_is_read_by_dna_variant_records_loader() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let dna_variant_records: Vec<DNAVariantRecord> = vec![
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 2,
            chromosome_1: "chr17".into(),
            position_1: 7674224,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674226,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read_1;read_2".into(),
            num_consensus_read_names: 2,
            read_names: "read_1;read_2".into(),
            num_read_names: 2
        },
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 1,
            chromosome_1: "chr17".into(),
            position_1: 7674200,
            strand_1: "-".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674231,
            strand_2: "-".into(),
            operation_2: "U".into(),
            sequence: "".into(),
            variant_size: Some(30),
            variant_type: "DEL".into(),
            consensus_read_names: "read_3".into(),
            num_consensus_read_names: 1,
            read_names: "read_3;read_4".into(),
            num_read_names: 2
        }
    ];
    let mut df_variant_records = dna_variant_records_to_dataframe(dna_variant_records.clone());

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    let mut tsv = fs::File::create(file.path()).unwrap();
    CsvWriter::new(&mut tsv)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df_variant_records)
        .unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let mut df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    let mut annotated_tsv = fs::File::create(annotated_file.path()).unwrap();
    CsvWriter::new(&mut annotated_tsv)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df_annotated)
        .unwrap();

    // The loader of integrate-vars reads the annotated file and returns the records that went in
    let loaded_dna_variant_records: Vec<DNAVariantRecord> = load_dna_variant_records(annotated_file.path().to_str().unwrap());
    assert_eq!(loaded_dna_variant_records, dna_variant_records);

    // The annotation set reads the same file
    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());
    assert_eq!(variant_call_annotation_set.annotations.len(), 2);
    assert_eq!(variant_call_annotation_set.get(1).variant_type, VariantType::Deletion);
    assert_eq!(&*variant_call_annotation_set.get(1).variant_sequence, "");
    assert_eq!(variant_call_annotation_set.get(1).position_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.get(2).position_2_annotation.reference_gene_ids.contains("ENSG00000141510.18"), true);
}

#[test]
fn append_variant_call_annotations_replaces_annotation_columns_of_annotated_table() {
    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-006-tumor_exacto_germline_dna_variants.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let df_annotated_1: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
    let df_annotated_2: DataFrame = append_variant_call_annotations(
        &df_annotated_1,
        &gene_annotator,
        1
    );

    assert_eq!(df_annotated_1.width(), df_variant_calls.width() + 8);
    assert_eq!(df_annotated_2.equals_missing(&df_annotated_1), true);
}

/// TP53 exon ENSE00003625790.1 is chr17:7675994-7676272. With protein-coding transcripts of
/// levels 1 and 2, the 20 bases on each side of it are intronic.
#[test]
fn append_variant_call_annotations_annotates_bases_between_positions() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    // variant_id 1: SNV on the first exonic base
    // variant_id 2: SNV on the last intronic base before the exon
    // variant_id 3: SNV on the last exonic base
    // variant_id 4: SNV on the first intronic base after the exon
    // variant_id 5: deletion of the first exonic base
    // variant_id 6: deletion of the last intronic base before the exon
    // variant_id 7: MNV on the last two exonic bases and the first intronic base
    // variant_id 8: deletion of the whole exon
    // variant_id 9: insertion between the last intronic base and the first exonic base
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t7675993\tchr17\t7675995\tSNV\tA\n\
         2\tchr17\t7675992\tchr17\t7675994\tSNV\tA\n\
         3\tchr17\t7676271\tchr17\t7676273\tSNV\tA\n\
         4\tchr17\t7676272\tchr17\t7676274\tSNV\tA\n\
         5\tchr17\t7675993\tchr17\t7675995\tDEL\t\n\
         6\tchr17\t7675992\tchr17\t7675994\tDEL\t\n\
         7\tchr17\t7676270\tchr17\t7676274\tMNV\tACG\n\
         8\tchr17\t7675974\tchr17\t7676292\tDEL\t\n\
         9\tchr17\t7675993\tchr17\t7675994\tINS\tACGT\n"
    ).unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1u8, 2u8])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1u8, 2u8]))
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    let mut genic_regions: Vec<[String; 4]> = Vec::new();
    for i in 0..df_annotated.height() {
        genic_regions.push([
            df_annotated.column("position_1_genic_region").unwrap().str().unwrap().get(i).unwrap().to_string(),
            df_annotated.column("position_1_plus_1_genic_region").unwrap().str().unwrap().get(i).unwrap().to_string(),
            df_annotated.column("position_2_minus_1_genic_region").unwrap().str().unwrap().get(i).unwrap().to_string(),
            df_annotated.column("position_2_genic_region").unwrap().str().unwrap().get(i).unwrap().to_string()
        ]);
    }

    //                            position_1  position_1 + 1  position_2 - 1  position_2
    assert_eq!(genic_regions[0], ["intronic", "exonic",       "exonic",       "exonic"]);
    assert_eq!(genic_regions[1], ["intronic", "intronic",     "intronic",     "exonic"]);
    assert_eq!(genic_regions[2], ["exonic",   "exonic",       "exonic",       "intronic"]);
    assert_eq!(genic_regions[3], ["exonic",   "intronic",     "intronic",     "intronic"]);
    assert_eq!(genic_regions[4], ["intronic", "exonic",       "exonic",       "exonic"]);
    assert_eq!(genic_regions[5], ["intronic", "intronic",     "intronic",     "exonic"]);
    assert_eq!(genic_regions[6], ["exonic",   "exonic",       "intronic",     "intronic"]);
    assert_eq!(genic_regions[7], ["intronic", "intronic",     "intronic",     "intronic"]);
    assert_eq!(genic_regions[8], ["intronic", "exonic",       "intronic",     "exonic"]);

    // The exon is named where the base lies in it, and only there
    let col_position_1_annotation = df_annotated.column("position_1_annotation").unwrap().str().unwrap();
    let col_position_1_plus_1_annotation = df_annotated.column("position_1_plus_1_annotation").unwrap().str().unwrap();
    let col_position_2_minus_1_annotation = df_annotated.column("position_2_minus_1_annotation").unwrap().str().unwrap();
    assert_eq!(col_position_1_annotation.get(0).unwrap().contains("ENST00000269305.9-ENSE00003625790.1"), false);
    assert_eq!(col_position_1_plus_1_annotation.get(0).unwrap().contains("ENST00000269305.9-ENSE00003625790.1"), true);
    assert_eq!(col_position_2_minus_1_annotation.get(0).unwrap().contains("ENST00000269305.9-ENSE00003625790.1"), true);
    assert_eq!(col_position_1_plus_1_annotation.get(1).unwrap().contains("ENST00000269305.9-ENSE00003625790.1"), false);
}

#[test]
fn append_variant_call_annotations_annotates_position_zero() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t0\tchr17\t0\tBND\tA\n"
    ).unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(df_annotated.column("position_1_plus_1_genic_region").unwrap().str().unwrap().get(0).unwrap(), "intergenic");
    assert_eq!(df_annotated.column("position_2_minus_1_genic_region").unwrap().str().unwrap().get(0).unwrap(), "intergenic");
}

#[test]
#[should_panic(expected = "variant_id 1 is on more than one row of df_variant_calls (row 3 is the second).")]
fn append_variant_call_annotations_panics_when_variant_id_is_on_two_rows() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t7674224\tchr17\t7674226\tSNV\tA\n\
         2\tchr17\t7675993\tchr17\t7675995\tSNV\tA\n\
         1\tchr18\t5170100\tchr18\t5170102\tSNV\tC\n"
    ).unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

#[test]
#[should_panic(expected = "variant_id 1 is on more than one row of df_variant_calls (row 3 is the second).")]
fn annotate_variant_calls_panics_when_variant_id_is_on_two_rows() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t7674224\tchr17\t7674226\tSNV\tA\n\
         2\tchr17\t7675993\tchr17\t7675995\tSNV\tA\n\
         1\tchr18\t5170100\tchr18\t5170102\tSNV\tC\n"
    ).unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
        .unwrap()
        .finish()
        .unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    annotate_variant_calls(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

#[test]
fn annotate_variant_calls_keeps_order_of_rows() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         7\tchr17\t7674224\tchr17\t7674226\tSNV\tA\n\
         3\tchr17\t7675993\tchr17\t7675995\tSNV\tA\n\
         9\tchr17\t7676271\tchr17\t7676273\tSNV\tA\n\
         1\tchr17\t7676272\tchr17\t7676274\tSNV\tA\n\
         8\tchr17\t7675974\tchr17\t7676292\tDEL\t\n\
         2\tchr17\t7675993\tchr17\t7675994\tINS\tACGT\n\
         10\tchr18\t5170100\tchr18\t5170102\tSNV\tC\n\
         5\tchr17\t7676155\tchr18\t5170101\tTRA\tG\n"
    ).unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variant_calls = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(file.path().to_str().unwrap().into()))
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
        4
    );

    let variant_ids: Vec<u64> = variant_call_annotation_set
        .to_dataframe()
        .column("variant_id").unwrap()
        .u64().unwrap()
        .into_no_null_iter()
        .collect();

    assert_eq!(variant_ids, vec![7, 3, 9, 1, 8, 2, 10, 5]);
}

#[test]
fn append_variant_call_annotations_appends_columns_to_table_without_rows() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    // A header and no row
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\torigin\tchromosome_1\tposition_1\tstrand_1\toperation_1\tchromosome_2\tposition_2\tstrand_2\toperation_2\tsequence\tvariant_size\tvariant_type\tconsensus_read_names\tnum_consensus_read_names\tread_names\tnum_read_names\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    let column_names: Vec<String> = df_annotated
        .get_column_names()
        .iter()
        .map(|name| name.to_string())
        .collect();
    assert_eq!(df_variant_calls.shape(), (0, 17));
    assert_eq!(df_annotated.shape(), (0, 25));
    assert_eq!(
        column_names,
        vec![
            "variant_id", "origin", "chromosome_1", "position_1", "strand_1", "operation_1",
            "chromosome_2", "position_2", "strand_2", "operation_2", "sequence", "variant_size",
            "variant_type", "consensus_read_names", "num_consensus_read_names", "read_names",
            "num_read_names",
            "position_1_genic_region", "position_1_annotation",
            "position_1_plus_1_genic_region", "position_1_plus_1_annotation",
            "position_2_minus_1_genic_region", "position_2_minus_1_annotation",
            "position_2_genic_region", "position_2_annotation"
        ]
    );
}

/// The variant caller writes a file of zero bytes when it finds no variant.
#[test]
fn append_variant_call_annotations_appends_columns_to_table_of_file_of_zero_bytes() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    let dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    write_tsv_file(dna_variant_records, file.path()).unwrap();
    assert_eq!(fs::metadata(file.path()).unwrap().len(), 0);

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let mut df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    let column_names: Vec<String> = df_annotated
        .get_column_names()
        .iter()
        .map(|name| name.to_string())
        .collect();
    assert_eq!(df_variant_calls.shape(), (0, 0));
    assert_eq!(df_annotated.shape(), (0, 8));
    assert_eq!(
        column_names,
        vec![
            "position_1_genic_region", "position_1_annotation",
            "position_1_plus_1_genic_region", "position_1_plus_1_annotation",
            "position_2_minus_1_genic_region", "position_2_minus_1_annotation",
            "position_2_genic_region", "position_2_annotation"
        ]
    );

    // The annotated table is written with its header
    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    let mut annotated_tsv = fs::File::create(annotated_file.path()).unwrap();
    CsvWriter::new(&mut annotated_tsv)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df_annotated)
        .unwrap();
    assert_eq!(
        fs::read_to_string(annotated_file.path()).unwrap(),
        "position_1_genic_region\tposition_1_annotation\tposition_1_plus_1_genic_region\tposition_1_plus_1_annotation\tposition_2_minus_1_genic_region\tposition_2_minus_1_annotation\tposition_2_genic_region\tposition_2_annotation\n"
    );
}

#[test]
fn annotate_variant_calls_returns_empty_set_for_table_without_rows() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

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

    assert_eq!(variant_call_annotation_set.annotations.len(), 0);
    assert_eq!(variant_call_annotation_set.to_dataframe().shape(), (0, 15));
}

#[test]
#[should_panic(expected = "chromosome_1 '17' on row 2 of df_variant_calls does not start with 'chr'.")]
fn append_variant_call_annotations_panics_when_chromosome_1_does_not_start_with_chr() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t7674224\tchr17\t7674226\tSNV\tA\n\
         2\t17\t7675993\tchr17\t7675995\tSNV\tA\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

#[test]
#[should_panic(expected = "chromosome_2 '18' on row 1 of df_variant_calls does not start with 'chr'.")]
fn append_variant_call_annotations_panics_when_chromosome_2_does_not_start_with_chr() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr17\t7676155\t18\t5170101\tTRA\tG\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

/// Every name among the first 150 rows is a number and the name of the last row is not.
#[test]
#[should_panic(expected = "chromosome_1 '17' on row 1 of df_variant_calls does not start with 'chr'.")]
fn append_variant_call_annotations_panics_when_chromosome_names_are_numbers() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let mut tsv: String = "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n".to_string();
    for variant_id in 1..=150 {
        tsv.push_str(&format!("{}\t17\t{}\t17\t{}\tSNV\tA\n", variant_id, 7674224 + variant_id, 7674226 + variant_id));
    }
    tsv.push_str("151\tX\t100\tX\t102\tSNV\tA\n");
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(file.path(), tsv).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());
    assert_eq!(df_variant_calls.height(), 151);
    assert_eq!(df_variant_calls.column("chromosome_1").unwrap().dtype(), &DataType::String);

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

/// The table of the variant caller holds IDs and positions as unsigned integers.
#[test]
fn append_variant_call_annotations_reads_table_of_variant_caller_that_was_not_written_to_file() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let df_variant_calls: DataFrame = dna_variant_records_to_dataframe(vec![
        DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: 7,
            chromosome_1: "chr17".into(),
            position_1: 7674224,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7674226,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read_1;read_2".into(),
            num_consensus_read_names: 2,
            read_names: "read_1;read_2".into(),
            num_read_names: 2
        }
    ]);
    assert_eq!(df_variant_calls.column("variant_id").unwrap().dtype(), &DataType::UInt32);
    assert_eq!(df_variant_calls.column("position_1").unwrap().dtype(), &DataType::UInt32);

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    assert_eq!(df_annotated.shape(), (1, 25));
    assert_eq!(df_annotated.column("variant_id").unwrap().dtype(), &DataType::UInt32);
    assert_eq!(df_annotated.column("position_1_plus_1_genic_region").unwrap().str().unwrap().get(0).unwrap(), "exonic");
}

/// The gene annotation of the tests has chr17 and chr18.
#[test]
fn get_chromosomes_missing_from_gene_annotation_returns_missing_chromosomes() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let col_chromosome_1: StringChunked = StringChunked::new("chromosome_1".into(), &["chr2", "chr17", "chr17", "chr1"]);
    let col_chromosome_2: StringChunked = StringChunked::new("chromosome_2".into(), &["chr2", "chr17", "chr18", "chrX"]);

    let missing_chromosomes: Vec<Box<str>> = get_chromosomes_missing_from_gene_annotation(
        &col_chromosome_1,
        &col_chromosome_2,
        &gene_annotator
    );

    assert_eq!(missing_chromosomes, vec!["chr1".into(), "chr2".into(), "chrX".into()]);
}

#[test]
fn get_chromosomes_missing_from_gene_annotation_returns_none_when_gene_annotation_has_all() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let col_chromosome_1: StringChunked = StringChunked::new("chromosome_1".into(), &["chr17", "chr17"]);
    let col_chromosome_2: StringChunked = StringChunked::new("chromosome_2".into(), &["chr17", "chr18"]);

    let missing_chromosomes: Vec<Box<str>> = get_chromosomes_missing_from_gene_annotation(
        &col_chromosome_1,
        &col_chromosome_2,
        &gene_annotator
    );

    assert_eq!(missing_chromosomes.len(), 0);
}

#[test]
#[should_panic(expected = "df_variant_calls and the gene annotation share no chromosome name. df_variant_calls has chr1, chrX. The gene annotation has chr17, chr18.")]
fn append_variant_call_annotations_panics_when_no_chromosome_is_in_gene_annotation() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchrX\t7674224\tchrX\t7674226\tSNV\tA\n\
         2\tchr1\t7675993\tchr1\t7675995\tSNV\tA\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );
}

#[test]
fn append_variant_call_annotations_annotates_table_when_one_chromosome_is_not_in_gene_annotation() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\n\
         1\tchr1\t7674224\tchr1\t7674226\tSNV\tA\n\
         2\tchr17\t7674224\tchr17\t7674226\tSNV\tA\n"
    ).unwrap();

    let df_variant_calls: DataFrame = read_variant_calls_tsv_file(file.path().to_str().unwrap());

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let df_annotated: DataFrame = append_variant_call_annotations(
        &df_variant_calls,
        &gene_annotator,
        1
    );

    let col_position_1_genic_region = df_annotated.column("position_1_genic_region").unwrap().str().unwrap();
    let col_position_1_annotation = df_annotated.column("position_1_annotation").unwrap().str().unwrap();
    assert_eq!(df_annotated.height(), 2);
    assert_eq!(col_position_1_genic_region.get(0).unwrap(), "intergenic");
    assert_eq!(col_position_1_annotation.get(0).unwrap(), "||");
    assert_eq!(col_position_1_genic_region.get(1).unwrap(), "exonic");
}
