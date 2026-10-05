use exacto_core::prelude::{Gencode, GenicRegion};
use std::fs;
use std::path::Path;
use tempfile::NamedTempFile;

use super::*;


#[test]
fn variant_call_annotation_set_read_tsv_file_returns_matches_1() {
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

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    annotate_variant_calls(&df_variant_calls, &gene_annotator, 1).to_tsv_file(annotated_file.path().to_str().unwrap());

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set.annotations.keys().len(), 1);
    assert_eq!(variant_call_annotation_set.annotations.get(&1).unwrap().position_1_annotation.genic_region, GenicRegion::Exonic);
}

#[test]
fn variant_call_annotation_set_read_tsv_file_returns_matches_2() {
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

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    annotate_variant_calls(&df_variant_calls, &gene_annotator, 1).to_tsv_file(annotated_file.path().to_str().unwrap());

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set.get(1).id, 1);
    assert_eq!(variant_call_annotation_set.get_by_range("chr17".into(), 7_000_000, 8_000_000).len(), 1);
    assert_eq!(variant_call_annotation_set.get_by_reference_transcript("ENST00000269305.9".into()).is_some(), true);
}

#[test]
fn variant_call_annotation_set_clone_equals_original() {
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

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    annotate_variant_calls(&df_variant_calls, &gene_annotator, 1).to_tsv_file(annotated_file.path().to_str().unwrap());

    let variant_call_annotation_set_1: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());
    let variant_call_annotation_set_2: VariantCallAnnotationSet = variant_call_annotation_set_1.clone();

    assert_eq!(variant_call_annotation_set_1, variant_call_annotation_set_2);
}

#[test]
fn variant_call_annotation_set_to_tsv_file_returns_matches() {
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

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    annotate_variant_calls(&df_variant_calls, &gene_annotator, 1).to_tsv_file(annotated_file.path().to_str().unwrap());

    let variant_call_annotation_set_1: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    variant_call_annotation_set_1.to_tsv_file(file.path().to_str().unwrap());

    let variant_call_annotation_set_2: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set_1, variant_call_annotation_set_2);
}

#[test]
fn variant_call_annotation_set_read_tsv_file_reads_annotations_of_bases_between_positions() {
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

    let annotated_file: NamedTempFile = NamedTempFile::new().unwrap();
    annotate_variant_calls(&df_variant_calls, &gene_annotator, 1).to_tsv_file(annotated_file.path().to_str().unwrap());

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(annotated_file.path().to_str().unwrap());

    // The SNV alters chr17:7674225, which is position_1 + 1 and position_2 - 1
    assert_eq!(variant_call_annotation_set.get(1).position_1, 7674224);
    assert_eq!(variant_call_annotation_set.get(1).position_2, 7674226);
    assert_eq!(variant_call_annotation_set.get(1).position_1_plus_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.get(1).position_2_minus_1_annotation.genic_region, GenicRegion::Exonic);
    assert_eq!(variant_call_annotation_set.get(1).position_1_plus_1_annotation, variant_call_annotation_set.get(1).position_2_minus_1_annotation);
    assert_eq!(variant_call_annotation_set.get(1).position_1_plus_1_annotation.reference_exon_ids.get("ENST00000269305.9").unwrap().as_ref(), "ENSE00003712342.1");
}

#[test]
#[should_panic(expected = "variant call ID 1 is already in the variant call annotation set.")]
fn variant_call_annotation_set_add_annotation_panics_when_variant_call_id_is_in_set() {
    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();

    variant_call_annotation_set.add_annotation(VariantCallAnnotation::new(
        1,
        "chr17".into(),
        7674224,
        "chr17".into(),
        7674226,
        VariantType::SingleNucleotideVariant,
        "A".into(),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic)
    ));
    variant_call_annotation_set.add_annotation(VariantCallAnnotation::new(
        1,
        "chr18".into(),
        5170100,
        "chr18".into(),
        5170102,
        VariantType::SingleNucleotideVariant,
        "C".into(),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic),
        PositionAnnotation::new(GenicRegion::Intergenic)
    ));
}

#[test]
#[should_panic(expected = "variant call ID 1 is already in the variant call annotation set.")]
fn variant_call_annotation_set_read_tsv_file_panics_when_variant_id_is_on_two_rows() {
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\tposition_1_genic_region\tposition_1_annotation\tposition_1_plus_1_genic_region\tposition_1_plus_1_annotation\tposition_2_minus_1_genic_region\tposition_2_minus_1_annotation\tposition_2_genic_region\tposition_2_annotation\n\
         1\tchr17\t100\tchr17\t102\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         1\tchr18\t200\tchr18\t202\tSNV\tC\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n"
    ).unwrap();

    VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());
}

#[test]
fn variant_call_annotation_set_to_dataframe_returns_rows_in_order_of_addition() {
    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    for variant_call_id in [7, 3, 9, 1, 8, 2, 10, 5, 4, 6] {
        variant_call_annotation_set.add_annotation(VariantCallAnnotation::new(
            variant_call_id,
            "chr17".into(),
            1000 + (variant_call_id as u32),
            "chr17".into(),
            1002 + (variant_call_id as u32),
            VariantType::SingleNucleotideVariant,
            "A".into(),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic)
        ));
    }

    let variant_ids: Vec<u64> = variant_call_annotation_set
        .to_dataframe()
        .column("variant_id").unwrap()
        .u64().unwrap()
        .into_no_null_iter()
        .collect();

    assert_eq!(variant_ids, vec![7, 3, 9, 1, 8, 2, 10, 5, 4, 6]);
}

#[test]
fn variant_call_annotation_set_read_tsv_file_keeps_order_of_rows() {
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    fs::write(
        file.path(),
        "variant_id\tchromosome_1\tposition_1\tchromosome_2\tposition_2\tvariant_type\tsequence\tposition_1_genic_region\tposition_1_annotation\tposition_1_plus_1_genic_region\tposition_1_plus_1_annotation\tposition_2_minus_1_genic_region\tposition_2_minus_1_annotation\tposition_2_genic_region\tposition_2_annotation\n\
         7\tchr17\t100\tchr17\t102\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         3\tchr17\t200\tchr17\t202\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         9\tchr17\t300\tchr17\t302\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         1\tchr17\t400\tchr17\t402\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         8\tchr17\t500\tchr17\t502\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n\
         2\tchr17\t600\tchr17\t602\tSNV\tA\tintergenic\t||\tintergenic\t||\tintergenic\t||\tintergenic\t||\n"
    ).unwrap();

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());

    let variant_ids: Vec<u64> = variant_call_annotation_set
        .to_dataframe()
        .column("variant_id").unwrap()
        .u64().unwrap()
        .into_no_null_iter()
        .collect();

    assert_eq!(variant_ids, vec![7, 3, 9, 1, 8, 2]);
}

#[test]
fn variant_call_annotation_set_get_by_range_returns_in_order_of_addition() {
    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    for variant_call_id in [7, 3, 9, 1, 8, 2, 10, 5, 4, 6] {
        variant_call_annotation_set.add_annotation(VariantCallAnnotation::new(
            variant_call_id,
            "chr17".into(),
            1000 + (variant_call_id as u32),
            "chr17".into(),
            1002 + (variant_call_id as u32),
            VariantType::SingleNucleotideVariant,
            "A".into(),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic),
            PositionAnnotation::new(GenicRegion::Intergenic)
        ));
    }

    // Variant call 1 has position_2 = 1003 and variant call 10 has position_1 = 1010
    let variant_call_ids: Vec<usize> = variant_call_annotation_set
        .get_by_range("chr17".into(), 1003, 1010)
        .iter()
        .map(|variant_call_annotation| variant_call_annotation.id)
        .collect();

    assert_eq!(variant_call_ids, vec![7, 3, 9, 1, 8, 2, 10, 5, 4, 6]);
}

#[test]
fn variant_call_annotation_set_get_by_reference_transcript_returns_in_order_of_addition() {
    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    for variant_call_id in [7, 3, 9, 1, 8, 2, 10, 5, 4, 6] {
        let mut position_annotation: PositionAnnotation = PositionAnnotation::new(GenicRegion::Intronic);
        position_annotation.add_reference_gene_id("ENSG001".into());
        position_annotation.add_reference_transcript_id("ENSG001".into(), "ENST001".into());
        variant_call_annotation_set.add_annotation(VariantCallAnnotation::new(
            variant_call_id,
            "chr17".into(),
            1000 + (variant_call_id as u32),
            "chr17".into(),
            1002 + (variant_call_id as u32),
            VariantType::SingleNucleotideVariant,
            "A".into(),
            position_annotation.clone(),
            position_annotation.clone(),
            position_annotation.clone(),
            position_annotation.clone()
        ));
    }

    let variant_call_ids: Vec<usize> = variant_call_annotation_set
        .get_by_reference_transcript("ENST001".into()).unwrap()
        .iter()
        .map(|variant_call_annotation| variant_call_annotation.id)
        .collect();

    assert_eq!(variant_call_ids, vec![7, 3, 9, 1, 8, 2, 10, 5, 4, 6]);
}

#[test]
fn variant_call_annotation_set_read_tsv_file_returns_empty_set_for_file_without_rows() {
    let variant_call_annotation_set_1: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    variant_call_annotation_set_1.to_tsv_file(file.path().to_str().unwrap());
    assert_eq!(fs::read_to_string(file.path()).unwrap().lines().count(), 1);

    let variant_call_annotation_set_2: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set_2.annotations.len(), 0);
    assert_eq!(variant_call_annotation_set_2, variant_call_annotation_set_1);
}

#[test]
fn variant_call_annotation_set_read_tsv_file_returns_empty_set_for_file_of_zero_bytes() {
    let file: NamedTempFile = NamedTempFile::new().unwrap();
    assert_eq!(fs::metadata(file.path()).unwrap().len(), 0);

    let variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::read_tsv_file(file.path().to_str().unwrap());

    assert_eq!(variant_call_annotation_set.annotations.len(), 0);
}
