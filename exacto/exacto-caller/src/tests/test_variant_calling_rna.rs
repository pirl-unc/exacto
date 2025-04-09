extern crate csv;
extern crate noodles_bam;
extern crate polars;

use bimap::BiMap;
use csv::ReaderBuilder;
use exacto_util::prelude::*;
use noodles_bam as bam;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_identify_closest_reference_transcript_id_1() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);

    let records_map = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_104905_566269/1/ccs"; // normal (unaltered) ENST00000269305.9 read
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let records: &Vec<bam::Record> = records_map.get(&read_id).unwrap();
    let read_sequence: Box<str> =  get_original_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_original_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());

    let alignment: Alignment = Alignment::new(
        read_id,
        read_sequence,
        quality_scores,
        records.clone()
    );

    // Identify exons
    let exons: Vec<TranscriptModelExon> = alignment.identify_exons(25);

    // Identify splice junctions
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = alignment.identify_splice_junctions(
        &chromosome_names_map,
        reference_genome_fasta_file,
        25
    );

    // Identify overlapping gene IDs
    let mut reference_gene_ids: HashSet<Box<str>> = Alignment::identify_overlapping_gene_ids(
        &exons,
        &gene_annotator,
        &chromosome_names_map
    );

    // Identify closest reference transcript IDs
    let reference_transcript_ids: Vec<Box<str>> = Alignment::identify_closest_reference_transcript_ids(
        &exons,
        &gene_annotator,
        &chromosome_names_map,
        &reference_gene_ids
    );

    assert!(reference_transcript_ids.len() == 1);
    assert!(reference_transcript_ids[0] == "ENST00000269305.9".into());
}

#[test]
fn test_identify_rna_variants_1() {
    let bam_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-100-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-100-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_2() {
    let bam_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-101-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-101-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}


#[test]
fn test_identify_rna_variants_3() {
    let bam_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-102-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-102-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_4() {
    let bam_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-103-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 2);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-103-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_5() {
    let bam_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-104-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-104-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_6() {
    let bam_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-105-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-105-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_7() {
    let bam_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-106-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-106-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_8() {
    let bam_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-107-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 2);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-107-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}

#[test]
fn test_identify_rna_variants_9() {
    let bam_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-108-tumor_minimap2_mdtagged_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new("src/tests/data/gtf/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38"
    );

    let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        25,
        25f32,
        1
    );

    assert!(transcript_model_set.transcript_models.iter().len() == 1);
    assert!(transcript_model_set.transcript_models.iter().next().unwrap().variant_calls.len() == 1);

    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    transcript_model_set.load_chromosome_names(chromosome_names_map);
    transcript_model_set.load_read_names(read_names_map);

    // Compare against the ground truth
    let (df_ref_matches,df_exons,df_splice_junctions,df_variant_calls) = transcript_model_set.to_dataframes(1);
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-108-tumor_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let chromosome_1: String = record[1].to_string();
        let position_1: i64 = record[2].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_1: String = record[4].to_string();
        let chromosome_2: String = record[5].to_string();
        let position_2: i64 = record[6].parse::<i64>().expect("Failed to convert &str to usize");
        let operation_2: String = record[8].to_string();
        let variant_type: String = record[10].to_string();
        let variant_sequence: String = record[11].to_string();

        let df_variants_: PolarsResult<DataFrame>;

        df_variants_ = df_variant_calls
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_1.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_2.to_string())))
            .filter(col("operation_1").eq(lit(operation_1.to_string())))
            .filter(col("operation_2").eq(lit(operation_2.to_string())))
            .filter(col("position_1").gt_eq(lit(position_1 - 100)))
            .filter(col("position_1").lt_eq(lit(position_1 + 100)))
            .filter(col("position_2").gt_eq(lit(position_2 - 100)))
            .filter(col("position_2").lt_eq(lit(position_2 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();

        assert!(df_variants_.unwrap().height() == 1);
    }
}