use exacto_core::prelude::*;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_genome_variation_graph_1() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_1.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATACGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_2() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_2.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new(),
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATACGTTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_3() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_3.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCACGTACAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_4() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_4.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTTTCC".into());
}

#[test]
fn test_genome_variation_graph_5() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_5.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTACGTAGGTAGCTAGCCGTACGTAGGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_6() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_6.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTCGTACGTAGCTAGCTACTAG".into());
}

#[test]
fn test_genome_variation_graph_7() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_7.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGCACCGCACGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_8() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_8.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "TAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_9() {
    let fasta_path = Path::new("src/tests/data/fasta/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_9.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "TACGCGTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_10() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_10.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 3);

    let mut sequence_exists_1: bool = false;
    let mut sequence_exists_2: bool = false;
    let mut sequence_exists_3: bool = false;

    for path in paths.iter() {
        if path.get_sequence() == "GCTCG".into() {
            sequence_exists_1 = true;
        }
        if path.get_sequence() == "TTTGC".into() {
            sequence_exists_2 = true;
        }
        if path.get_sequence() == "ACGCGTACGTGTACGTAGCTACCCTTTGGGAAACGC".into() ||
            reverse_complement(&*path.get_sequence()) == "ACGCGTACGTGTACGTAGCTACCCTTTGGGAAACGC".into() {
            sequence_exists_3 = true;
        }
    }

    assert!(sequence_exists_1);
    assert!(sequence_exists_2);
    assert!(sequence_exists_3);
}

#[test]
fn test_genome_variation_graph_11() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_11.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_node_ids(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 4);

    let mut sequence_exists_1: bool = false;
    let mut sequence_exists_2: bool = false;
    let mut sequence_exists_3: bool = false;
    let mut sequence_exists_4: bool = false;

    for path in paths.iter() {
        if path.get_sequence() == "ATGCGTACGTAGCTA".into() {
            sequence_exists_1 = true;
        }
        if path.get_sequence() == "GGGTTTCCCAAAGGG".into() {
            sequence_exists_2 = true;
        }
        if path.get_sequence() == "GGAAAGCTAG".into() ||
            reverse_complement(&*path.get_sequence()) == "GGAAAGCTAG".into() {
            sequence_exists_3 = true;
        }
        if path.get_sequence() == "GGATCGTATCTGACGTATGA".into() {
            sequence_exists_4 = true;
        }
    }

    assert!(sequence_exists_1);
    assert!(sequence_exists_2);
    assert!(sequence_exists_3);
    assert!(sequence_exists_4);
}

#[test]
fn test_genome_variation_graph_12() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_12.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGGATACGTAGCTAGGACTAG".into());
}

#[test]
fn test_genome_variation_graph_13() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_13.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCAGGTACGTAGCTAAGGCTAG".into());
}

#[test]
fn test_genome_variation_graph_14() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_14.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCTTAAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_15() {
    let fasta_path = Path::new("src/tests/data/fasta/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/sample_dna_variant_callset_15.tsv");
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

    let vargraphs: Vec<VarGraph> = build_genome_variation_graph(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        2
    );

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: HashSet<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATTCGTAGCTAGCTAG".into());
}

#[test]
fn test_transcriptome_variation_graph_1() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_1.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGATAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_2() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_2.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCCCCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_3() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_3.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGATAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_4() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_4.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTAGCAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_5() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_5.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCGTCGAGCACCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_6() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_6.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_7() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_7.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTCTCGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_8() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_8.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTCTCGCCGGTCGA".into());
}

#[test]
fn test_transcriptome_variation_graph_9() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_9.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAATCGC".into());
}

#[test]
fn test_transcriptome_variation_graph_10() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_10.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "TAACTGCGATTACTG".into());
}

#[test]
fn test_transcriptome_variation_graph_11() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_11.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "AGCTAAGCTA".into());
}

#[test]
fn test_transcriptome_variation_graph_12() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_12.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_13() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_13.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "TCGAAGCGTCGAGC".into());
}

#[test]
fn test_transcriptome_variation_graph_14() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_14.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGGGGGAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_15() {
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_15.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGCGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_16() {
    // Forward strand, 1-exon transcript (chrA:1-10, +)
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_16.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGTACGT".into());
}

#[test]
fn test_transcriptome_variation_graph_17() {
    // Reverse strand, 1-exon transcript (chrA:1-10, -)
    // Expected: RC("ATGCGTACGT") = "ACGTACGCAT"
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_17.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ACGTACGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_18() {
    // Forward strand, 3-exon transcript (chrA:1-5, 11-15, 21-25, all +)
    // Expected: "ATGCG" + "AGCTA" + "AGCGT" = "ATGCGAGCTAAGCGT"
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_18.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_19() {
    // Reverse strand, 3-exon transcript (chrA:21-25, 11-15, 1-5, all -)
    // Expected: RC("AGCGT") + RC("AGCTA") + RC("ATGCG") = "ACGCTTAGCTCGCAT"
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_19.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "ACGCTTAGCTCGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_20() {
    // Reverse strand, 2-exon transcript (chrA:11-15, 1-5, all -)
    // Expected: RC("AGCTA") + RC("ATGCG") = "TAGCTCGCAT"
    let fasta_path = Path::new("src/tests/data/fasta/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new("src/tests/data/tsv/sample_rna_variant_callset_20.tsv");
    let variant_tsv_full_path = fs::canonicalize(variant_tsv_path).unwrap();
    let variant_tsv_file: &str = variant_tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variant_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraphs: Vec<(usize, Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    );

    assert_eq!(vargraphs.len(), 1);

    let transcript_model_id: usize = vargraphs.get(0).unwrap().0;
    let reference_transcript_ids: &str = &*vargraphs.get(0).unwrap().1;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().2.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_model_id,
        reference_transcript_ids,
        &vargraph.get_variant_node_ids().into_iter().collect()
    );

    assert_eq!(path.get_sequence(), "TAGCTCGCAT".into());
}
