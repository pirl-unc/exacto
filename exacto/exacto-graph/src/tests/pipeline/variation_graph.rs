use std::collections::HashSet;
use std::fs;
use std::path::Path;
use super::*;


#[test]
fn test_genome_variation_graph_1() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_1.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATACGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_2() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_2.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new(),
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATACGTTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_3() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_3.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCACGTACAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_4() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_4.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTTTCC".into());
}

#[test]
fn test_genome_variation_graph_5() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_5.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTACGTAGGTAGCTAGCCGTACGTAGGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_6() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_6.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGTCGTACGTAGCTAGCTACTAG".into());
}

#[test]
fn test_genome_variation_graph_7() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_7.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGCACCGCACGTAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_8() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_8.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "TAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_9() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample2.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_9.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "TACGCGTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_10() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_10.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
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
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_11.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_node_ids(),
        &HashSet::new()
    );

    // chrC has no variant, so it has no graph; its reference is written by
    // find_genome_variation_graph_sequences instead.
    assert_eq!(paths.len(), 3);

    let mut sequence_exists_1: bool = false;
    let mut sequence_exists_2: bool = false;
    let mut sequence_exists_3: bool = false;

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
    }

    assert!(sequence_exists_1);
    assert!(sequence_exists_2);
    assert!(sequence_exists_3);
}

#[test]
fn test_genome_variation_graph_12() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_12.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCGGATACGTAGCTAGGACTAG".into());
}

#[test]
fn test_genome_variation_graph_13() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_13.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCAGGTACGTAGCTAAGGCTAG".into());
}

#[test]
fn test_genome_variation_graph_14() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_14.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCTTAAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_15() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_15.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    assert_eq!(paths.len(), 1);
    assert_eq!(paths.iter().next().unwrap().get_sequence(), "ATGCATTCGTAGCTAGCTAG".into());
}

#[test]
fn test_transcriptome_variation_graph_1() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_1.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGATAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_2() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_2.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCCCCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_3() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_3.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGATAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_4() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_4.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTAGCAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_5() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_5.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCGTCGAGCACCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_6() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_6.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_7() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_7.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTCTCGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_8() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_8.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTCTCGCCGGTCGA".into());
}

#[test]
fn test_transcriptome_variation_graph_9() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_9.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAATCGC".into());
}

#[test]
fn test_transcriptome_variation_graph_10() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_10.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "TAACTGCGATTACTG".into());
}

#[test]
fn test_transcriptome_variation_graph_11() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_11.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "AGCTAAGCTA".into());
}

#[test]
fn test_transcriptome_variation_graph_12() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_12.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_13() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_13.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "TCGAAGCGTCGAGC".into());
}

#[test]
fn test_transcriptome_variation_graph_14() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_14.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGGGGGAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_15() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_15.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGCGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_16() {
    // Forward strand, 1-exon transcript (chrA:1-10, +)
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_16.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGTACGT".into());
}

#[test]
fn test_transcriptome_variation_graph_17() {
    // Reverse strand, 1-exon transcript (chrA:1-10, -)
    // Expected: RC("ATGCGTACGT") = "ACGTACGCAT"
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_17.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ACGTACGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_18() {
    // Forward strand, 3-exon transcript (chrA:1-5, 11-15, 21-25, all +)
    // Expected: "ATGCG" + "AGCTA" + "AGCGT" = "ATGCGAGCTAAGCGT"
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_18.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ATGCGAGCTAAGCGT".into());
}

#[test]
fn test_transcriptome_variation_graph_19() {
    // Reverse strand, 3-exon transcript (chrA:21-25, 11-15, 1-5, all -)
    // Expected: RC("AGCGT") + RC("AGCTA") + RC("ATGCG") = "ACGCTTAGCTCGCAT"
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_19.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "ACGCTTAGCTCGCAT".into());
}

#[test]
fn test_transcriptome_variation_graph_20() {
    // Reverse strand, 2-exon transcript (chrA:11-15, 1-5, all -)
    // Expected: RC("AGCTA") + RC("ATGCG") = "TAGCTCGCAT"
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_20.tsv");
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

    let vargraphs: Vec<(Box<str>, VarGraph)> = build_transcriptome_variation_graph(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2
    ).unwrap();

    assert_eq!(vargraphs.len(), 1);

    let transcript_name: &str = &*vargraphs.get(0).unwrap().0;
    let vargraph: VarGraph = vargraphs.get(0).unwrap().1.clone();

    let path: VarGraphPath = vargraph.find_transcript_path(
        transcript_name,
        &vargraph.get_variant_node_ids().into_iter().collect()
    ).unwrap();

    assert_eq!(path.get_sequence(), "TAGCTCGCAT".into());
}

#[test]
fn test_genome_variation_graph_minus_strand_sequences_are_genome_forward() {
    // Review F-01: a variant's sequence is stored genome-forward, so strands "-" and "*" must
    // not reverse-complement it. The table also has no variant_type or num_cycles (F-02).
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_dna_variant_callset_16.tsv");
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
    ).unwrap();

    let vargraph: VarGraph = VarGraph::merge(vargraphs);

    let paths: Vec<VarGraphPath> = vargraph.find_genome_paths(
        &vargraph.get_variant_node_ids().into_iter().collect(),
        &HashSet::new()
    );

    // SNV G5A on "-" and an insertion of ACG after base 10 on "*"
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].get_sequence(), "ATGCATACGTACGAGCTAGCTAG".into());
}

#[test]
fn test_genome_variation_graph_reads_caller_insertion_on_minus_strand() {
    // Review F-01 and F-02: dna-002 as exacto-caller wrote it (no num_cycles column). Its
    // insertion is on "-" and is spelled genome-forward, as in the simulation truth.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/dna-002-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
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

    let mut sequences: Vec<(Box<str>, bool)> = Vec::new();
    find_genome_variation_graph_sequences(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        true,
        2,
        |sequence, is_variant| {
            sequences.push((sequence, is_variant));
            Ok(())
        }
    ).unwrap();

    // One chr17 with the insertion, then the reference of chr17 and chr18
    assert_eq!(sequences.len(), 3);
    assert_eq!(sequences.iter().map(|(_, is_variant)| *is_variant).collect::<Vec<bool>>(), vec![true, false, false]);
    let expected_junction: String = format!(
        "{}{}{}",
        get_fasta_sequence("chr17", 7674206, 7674225, fasta_file).to_uppercase(),
        "ACGTACGTGGTATGCATGCTGAGACTGAGG",
        get_fasta_sequence("chr17", 7674226, 7674245, fasta_file).to_uppercase()
    );
    assert!(sequences[0].0.to_uppercase().contains(&expected_junction));
    assert_eq!(sequences[0].0.len(), 10_000_000 + 30);
    assert_eq!(sequences[1].0, get_fasta_sequence("chr17", 1, 10_000_000, fasta_file));
    assert_eq!(sequences[2].0, get_fasta_sequence("chr18", 1, 10_000_000, fasta_file));
}

#[test]
fn test_genome_variation_graph_reads_caller_tandem_duplication() {
    // Review F-02: dna-008 as exacto-caller wrote it. The duplication is a U/D row with no
    // num_cycles, taken as 2 copies; the 1-base insertion T is on "-".
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/dna-008-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
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

    let mut sequences: Vec<Box<str>> = Vec::new();
    find_genome_variation_graph_sequences(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        false,
        2,
        |sequence, _| {
            sequences.push(sequence);
            Ok(())
        }
    ).unwrap();

    // chr17 up to the end of the duplicated span, the junction base G, the span again, the rest
    let expected: String = format!(
        "{}T{}G{}{}",
        get_fasta_sequence("chr17", 1, 7668420, fasta_file),
        get_fasta_sequence("chr17", 7668421, 7687489, fasta_file),
        get_fasta_sequence("chr17", 7668421, 7687489, fasta_file),
        get_fasta_sequence("chr17", 7687490, 10_000_000, fasta_file)
    );
    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].to_uppercase(), expected.to_uppercase());
}

#[test]
fn test_transcriptome_variation_graph_reads_caller_model_alignments() {
    // Review F-01 and F-03: TP53 (rna-100, minus strand) from the transcript model alignments
    // table of call-rna-transcript-vars. Every row's sequence is genome-forward, so the transcript
    // is the reverse complement of each row, taken in index order. Base 879 is the mismatch row.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/rna-100-tumor_minimap2_mdtagged_sorted_exacto_assembled_transcript_model_alignments.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let tsv_file: &str = tsv_full_path.to_str().unwrap();

    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_transcript_structures: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();
    find_transcriptome_variation_graph_sequences(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2,
        0,
        |name, sequence| {
            sequences.push((name, sequence));
            Ok(())
        }
    ).unwrap();

    let col_index = df_transcript_structures.column("index").unwrap().i64().unwrap();
    let col_sequence = df_transcript_structures.column("sequence").unwrap().str().unwrap();
    let mut rows: Vec<(i64, &str)> = (0..df_transcript_structures.height())
        .map(|i| (col_index.get(i).unwrap(), col_sequence.get(i).unwrap_or("")))
        .collect();
    rows.sort();
    let read: String = rows.iter().map(|(_, sequence)| reverse_complement(sequence).to_string()).collect();

    assert_eq!(sequences.len(), 1);
    assert_eq!(&*sequences[0].0, "m64012_507476_774164/1/ccs");
    assert_eq!(read.len(), 2512);
    assert_eq!(&read[878..879], "T");
    assert_eq!(sequences[0].1.to_uppercase(), read.to_uppercase());
}

#[test]
fn test_genome_variation_graph_sequences_cover_alleles_and_write_references() {
    // Review F-04: three sites with two alternative alleles each. Every allele is written once,
    // in 2 sequences rather than 2^3, and the reference of every contig follows.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let df_variants: DataFrame = polars::df!(
        "variant_id" => [1i64, 2, 3, 4, 5, 6],
        "chromosome_1" => ["chrA"; 6],
        "position_1" => [4i64, 4, 9, 9, 14, 14],
        "operation_1" => ["D"; 6],
        "strand_1" => ["+"; 6],
        "chromosome_2" => ["chrA"; 6],
        "position_2" => [6i64, 6, 11, 11, 16, 16],
        "operation_2" => ["U"; 6],
        "strand_2" => ["+"; 6],
        "sequence" => ["A", "C", "C", "A", "G", "T"]
    ).unwrap();

    let mut sequences: Vec<(Box<str>, bool)> = Vec::new();
    find_genome_variation_graph_sequences(
        fasta_file,
        &df_variants,
        VarGraphTypes::Population,
        true,
        2,
        |sequence, is_variant| {
            sequences.push((sequence, is_variant));
            Ok(())
        }
    ).unwrap();

    let variant_sequences: Vec<&str> = sequences.iter().filter(|(_, v)| *v).map(|(s, _)| &**s).collect();
    let reference_sequences: Vec<&str> = sequences.iter().filter(|(_, v)| !*v).map(|(s, _)| &**s).collect();
    assert_eq!(variant_sequences.len(), 2);
    for (position, alleles) in [(5usize, "AC"), (10, "AC"), (15, "GT")] {
        let mut seen: Vec<char> = variant_sequences.iter().map(|s| s.chars().nth(position - 1).unwrap()).collect();
        seen.sort();
        assert_eq!(seen.into_iter().collect::<String>(), alleles);
    }
    assert_eq!(reference_sequences, vec!["ATGCGTACGTAGCTAGCTAG", "GGGTTTCCCAAAGGGTTTCC", "GGATCGTATCTGACGTATGA"]);
}

#[test]
fn test_genome_variation_graph_sequences_follow_fasta_order() {
    // Review F-08: rows given chrC, chrA, chrB come out in FASTA order, the same in every run.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let df_variants: DataFrame = polars::df!(
        "variant_id" => [1i64, 2, 3],
        "chromosome_1" => ["chrC", "chrA", "chrB"],
        "position_1" => [4i64, 4, 4],
        "operation_1" => ["D"; 3],
        "strand_1" => ["+"; 3],
        "chromosome_2" => ["chrC", "chrA", "chrB"],
        "position_2" => [6i64, 6, 6],
        "operation_2" => ["U"; 3],
        "strand_2" => ["+"; 3],
        "sequence" => ["T", "A", "A"]
    ).unwrap();

    for _ in 0..5 {
        let mut sequences: Vec<Box<str>> = Vec::new();
        find_genome_variation_graph_sequences(
            fasta_file,
            &df_variants,
            VarGraphTypes::Individual,
            false,
            3,
            |sequence, _| {
                sequences.push(sequence);
                Ok(())
            }
        ).unwrap();
        assert_eq!(
            sequences,
            vec!["ATGCATACGTAGCTAGCTAG".into(), "GGGTATCCCAAAGGGTTTCC".into(), "GGATTGTATCTGACGTATGA".into()]
        );
    }
}

#[test]
fn test_transcriptome_variation_graph_skips_unresolved_model() {
    // Review F-09: model "good" has its rows out of index order; model "broken" lacks its first
    // junction row, so its rows make two paths; model "back_splice" has a U/D row and no
    // num_cycles. "good" is written and the other two are skipped.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample3.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let variant_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample_rna_variant_callset_21.tsv");
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

    let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();
    find_transcriptome_variation_graph_sequences(
        fasta_file,
        &df_transcript_structures,
        VarGraphTypes::Individual,
        2,
        1,
        |name, sequence| {
            sequences.push((name, sequence));
            Ok(())
        }
    ).unwrap();

    assert_eq!(sequences, vec![("good".into(), "ATGCGAGATAAGCGT".into())]);
}

#[test]
fn test_genome_variation_graph_does_not_stitch_across_chromosomes() {
    // Review F-10: a deletion chrA:5-15 and a translocation chrA:10 to chrB:6. chrB:6 equals the
    // deletion's position_1 + 1, but on another chromosome, so the two are not stitched.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let df_variants: DataFrame = polars::df!(
        "variant_id" => [1i64, 2],
        "chromosome_1" => ["chrA", "chrA"],
        "position_1" => [5i64, 10],
        "operation_1" => ["D", "D"],
        "strand_1" => ["+", "+"],
        "chromosome_2" => ["chrA", "chrB"],
        "position_2" => [15i64, 6],
        "operation_2" => ["U", "U"],
        "strand_2" => ["+", "+"],
        "sequence" => ["", ""]
    ).unwrap();

    let mut sequences: Vec<Box<str>> = Vec::new();
    find_genome_variation_graph_sequences(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        false,
        2,
        |sequence, _| {
            sequences.push(sequence);
            Ok(())
        }
    ).unwrap();
    sequences.sort();

    assert_eq!(sequences, vec!["ATGCGAGCTAG".into(), "TACGTTCCCAAAGGGTTTCC".into()]);
}

#[test]
fn test_genome_variation_graph_rejects_rows_outside_the_fasta() {
    // Review F-11: each bad row is an error that names the row, not a panic.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let cases: Vec<(&str, i64, &str, i64, &str, &str)> = vec![
        // chromosome, position_1, strand, position_2, expected error
        ("chr1", 4, "+", 6, "", "variant_id 1: chromosome_1 chr1 is not in the FASTA file"),
        ("chrA", 0, "+", 2, "", "variant_id 1: position_1 0 is not a base of chrA (1-20)"),
        ("chrA", 19, "+", 21, "", "variant_id 1: position_2 21 is not a base of chrA (1-20)"),
        ("chrA", 4, "x", 6, "", "variant_id 1: strand_1 'x' is not +, - or *"),
        ("chrA", 4, "+", 6, "2", "variant_id 1: num_cycles 2 is given for a row that is not a duplication")
    ];
    for (chromosome, position_1, strand, position_2, num_cycles, expected) in cases {
        let df_variants: DataFrame = polars::df!(
            "variant_id" => [1i64],
            "chromosome_1" => [chromosome],
            "position_1" => [position_1],
            "operation_1" => ["D"],
            "strand_1" => [strand],
            "chromosome_2" => [chromosome],
            "position_2" => [position_2],
            "operation_2" => ["U"],
            "strand_2" => ["+"],
            "sequence" => ["A"],
            "num_cycles" => [num_cycles]
        ).unwrap();
        let error = build_genome_variation_graph(fasta_file, &df_variants, VarGraphTypes::Individual, 1).unwrap_err();
        assert_eq!(error.to_string(), expected);
    }

    let df_variants: DataFrame = polars::df!(
        "variant_id" => [1i64],
        "chromosome_1" => ["chrA"],
        "position_1" => [4i64],
        "operation_1" => ["D"],
        "strand_1" => ["+"],
        "chromosome_2" => ["chrA"],
        "position_2" => [6i64],
        "operation_2" => ["U"],
        "strand_2" => ["+"]
    ).unwrap();
    let error = build_genome_variation_graph(fasta_file, &df_variants, VarGraphTypes::Individual, 1).unwrap_err();
    assert_eq!(error.to_string(), "column sequence: not in the table");
}

#[test]
fn test_genome_variation_graph_reads_numeric_chromosome_names() {
    // Review F-11: a chromosome named 17 is read by the CSV reader as a number.
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-graph/sample4.fa");
    let fasta_file_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_file_path.to_str().unwrap();

    let df_variants: DataFrame = polars::df!(
        "variant_id" => [1i64],
        "chromosome_1" => [17i64],
        "position_1" => [4i64],
        "operation_1" => ["D"],
        "strand_1" => ["+"],
        "chromosome_2" => [17i64],
        "position_2" => [6i64],
        "operation_2" => ["U"],
        "strand_2" => ["+"],
        "sequence" => ["A"]
    ).unwrap();

    let mut sequences: Vec<(Box<str>, bool)> = Vec::new();
    find_genome_variation_graph_sequences(
        fasta_file,
        &df_variants,
        VarGraphTypes::Individual,
        true,
        1,
        |sequence, is_variant| {
            sequences.push((sequence, is_variant));
            Ok(())
        }
    ).unwrap();

    assert_eq!(
        sequences,
        vec![
            ("ATGCATACGTAGCTAGCTAG".into(), true),
            ("ATGCGTACGTAGCTAGCTAG".into(), false),
            ("GGGTTTCCCAAAGGGTTTCC".into(), false)
        ]
    );
}
