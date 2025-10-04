use polars::prelude::*;
use std::fs;
use std::path::Path;
use crate::prelude::*;


#[test]
fn test_variation_graph_1() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_variant_callset_1.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variants: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraph: VarGraph = build_variation_graph(fasta_file, &df_variants);
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths.get(0).unwrap().get_sequence(), "ATGCATACGTAGCTAGCTAG".into());
}

#[test]
fn test_variation_graph_2() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_variant_callset_2.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variants: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraph: VarGraph = build_variation_graph(fasta_file, &df_variants);
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths.get(0).unwrap().get_sequence(), "ATGCATACGTTAGCTAG".into());
}

#[test]
fn test_variation_graph_3() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_variant_callset_3.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variants: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraph: VarGraph = build_variation_graph(fasta_file, &df_variants);
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths.get(0).unwrap().get_sequence(), "ATGCACGTACAGCTAGCTAG".into());
}

#[test]
fn test_variation_graph_4() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_variant_callset_4.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variants: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraph: VarGraph = build_variation_graph(fasta_file, &df_variants);
    let paths: Vec<VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());
    assert_eq!(paths.len(), 1);
    assert_eq!(paths.get(0).unwrap().get_sequence(), "ATGCGTTTCC".into());
}
