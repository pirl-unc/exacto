use std::collections::{HashMap, HashSet};
use csv::ReaderBuilder;
use exacto_core::prelude::*;
use polars::prelude::*;
use std::fs;
use std::fs::File;
use std::path::Path;
use tempfile::NamedTempFile;

use crate::prelude::*;


#[test]
fn test_variant_calling_rna_1() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();

    let file: NamedTempFile = NamedTempFile::new().unwrap();
    rna_variant_callset.to_tsv_file(file.path().to_str().unwrap(), 1);
    let rna_variant_callset_2: RNAVariantCallSet = RNAVariantCallSet::read_tsv_file(file.path().to_str().unwrap());
    assert_eq!(rna_variant_callset.to_dataframe(1), rna_variant_callset_2.to_dataframe(1));

    let variant_record: &VariantRecord = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1)
        .get(&vec!["ENST00000269305.9".into()])
        .unwrap()
        .first()
        .unwrap()
        .get_consensus_record()
        .0;
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(*variant_record.get_variant_type(), VariantType::SingleNucleotideVariant);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_2() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_record: &VariantRecord = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1)
        .get(&vec!["ENST00000269305.9".into()])
        .unwrap()
        .first()
        .unwrap()
        .get_consensus_record()
        .0;
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(*variant_record.get_variant_type(), VariantType::Insertion);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_3() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_record: &VariantRecord = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1)
        .get(&vec!["ENST00000269305.9".into()])
        .unwrap()
        .first()
        .unwrap()
        .get_consensus_record()
        .0;
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(*variant_record.get_variant_type(), VariantType::Deletion);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_4() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut fusion_gene_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000570791.5".into(), "ENST00000698746.1".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::FusionGene {
            fusion_gene_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(fusion_gene_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000570791.5,ENST00000698746.1");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_5() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut exon_truncation_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000269305.9".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::ExonTruncation {
            exon_truncation_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(exon_truncation_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_6() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut exon_truncation_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000269305.9".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::ExonTruncation {
            exon_truncation_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(exon_truncation_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_7() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut exon_truncation_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000269305.9".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::ExonTruncation {
            exon_truncation_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(exon_truncation_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_8() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut cryptic_exon_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000269305.9".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::CrypticExon {
            cryptic_exon_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(cryptic_exon_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_9() {
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut intron_retention_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000269305.9".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::IntronRetention {
            intron_retention_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(intron_retention_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_10() {
    let bam_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-109-tumor_minimap2_mdtagged_sorted.bam.bai");
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut fusion_gene_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000263087.9".into(), "ENST00000333813.4".into(), "ENST00000570791.5".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::FusionGene {
            fusion_gene_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(fusion_gene_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000263087.9,ENST00000333813.4,ENST00000570791.5");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-109-tumor_ground_truth.tsv");
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_11() {
    let bam_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-110-tumor_minimap2_mdtagged_sorted.bam.bai");
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut fusion_gene_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000250113.12".into(), "ENST00000263092.11".into(), "ENST00000355530.7".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::FusionGene {
            fusion_gene_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(fusion_gene_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000250113.12,ENST00000263092.11,ENST00000355530.7");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-110-tumor_ground_truth.tsv");
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
        df_variants_ = df_variant_calls_filtered
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
fn test_variant_calling_rna_12() {
    let bam_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new("src/tests/data/bam/rna-111-tumor_minimap2_mdtagged_sorted.bam.bai");
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
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
        bam_file,
        bam_bai_file,
        reference_genome_fasta_file,
        &gene_annotator,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.95f32,
        25,
        25,
        1
    );
    let rna_variant_callset: RNAVariantCallSet = transcript_model_set.get_variant_call_set();
    let variant_calls: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = rna_variant_callset
        .get_variant_calls_for_transcript_model_id(1);
    let mut circular_rna_found: bool = false;
    for variant_call in variant_calls.get(&vec!["ENST00000254719.10".into()]).unwrap().iter() {
        let variant_record: &VariantRecord = variant_call.get_consensus_record().0;
        if variant_record.get_variant_type() == &VariantType::CircularRNA {
            circular_rna_found = true;
        }
    }
    assert_eq!(rna_variant_callset.get_transcript_model_ids().len(), 1);
    assert_eq!(circular_rna_found, true);

    // Compare against the ground truth
    let mut df_variant_calls: DataFrame = rna_variant_callset.to_dataframe(1);
    let mask = df_variant_calls
        .column("reference_transcript_ids")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000254719.10");
    let df_variant_calls_filtered = df_variant_calls.filter(&mask).unwrap();
    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-111-tumor_ground_truth.tsv");
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
        let df_variants_1: PolarsResult<DataFrame> = df_variant_calls_filtered
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
        let df_variants_2: PolarsResult<DataFrame> = df_variant_calls_filtered
            .clone()
            .lazy()
            .filter(col("chromosome_1").eq(lit(chromosome_2.to_string())))
            .filter(col("chromosome_2").eq(lit(chromosome_1.to_string())))
            .filter(col("operation_1").eq(lit(operation_2.to_string())))
            .filter(col("operation_2").eq(lit(operation_1.to_string())))
            .filter(col("position_1").gt_eq(lit(position_2 - 100)))
            .filter(col("position_1").lt_eq(lit(position_2 + 100)))
            .filter(col("position_2").gt_eq(lit(position_1 - 100)))
            .filter(col("position_2").lt_eq(lit(position_1 + 100)))
            .filter(col("variant_type").lt_eq(lit(variant_type.to_string())))
            .collect();
        assert!(df_variants_1.unwrap().height() == 1 || df_variants_2.unwrap().height() == 1);
    }
}
