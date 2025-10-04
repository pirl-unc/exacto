use std::fs;
use std::path::Path;
use polars::prelude::*;
use crate::prelude::*;

// #[test]
// fn test_playground() {
//     let fasta_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/exacto/test/data/fasta/hg38_chr17_1-8M.fa";
//     let tsv_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_scga/mini/simulated_variants/scga-0001_germline_variants.tsv";
//
//     let parse_options = CsvParseOptions::default()
//         .with_separator(b'\t');
//     let df_variants: DataFrame = CsvReadOptions::default()
//         .with_parse_options(parse_options)
//         .with_has_header(true)
//         .try_into_reader_with_file_path(Some(tsv_file.into()))
//         .unwrap()
//         .finish()
//         .unwrap();
//
//     let vargraph: VarGraph = build_variation_graph(fasta_file, &df_variants);
// }
