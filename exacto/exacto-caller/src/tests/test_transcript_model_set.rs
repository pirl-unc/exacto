use csv::ReaderBuilder;
use exacto_core::prelude::*;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_transcript_model_set_1() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let mut df_transcript_structures_filtered: DataFrame = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-100-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_2() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-101-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_3() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-102-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        // Use >= 1 (not == 1): a BAM may contain multiple reads that all
        // align to the same reference transcript (e.g., two reads both
        // matching TP53). Each read becomes its own AssembledTranscript
        // and produces its own TranscriptModel for the same RT — so a
        // single ground-truth row can legitimately be present multiple
        // times in the dataframe, once per AT. transcript_model_id (the
        // distinguishing column) is non-deterministic across runs, so
        // we can't tighten this back to == 1 without redesigning the
        // ground-truth schema. The check still validates that the
        // ground-truth alignment row IS produced; it just doesn't
        // guarantee uniqueness.
        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_4() {
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
        1,
        1_000,
        ""
    );

    let structure_records: Vec<TranscriptModelStructureRecord> =
        build_transcript_model_structure_records(&transcript_model_set).collect();
    let expected_rt_set: HashSet<&str> = HashSet::from([
        "ENST00000698746.1",
        "ENST00000570791.5"
    ]);
    let records: Vec<&TranscriptModelStructureRecord> = structure_records
        .iter()
        .filter(|r| {
            let rt_set: HashSet<&str> = r.reference_transcript_id.split(',').collect();
            rt_set == expected_rt_set
        })
        .collect();
    let df_transcript_structures_filtered: DataFrame = transcript_model_structure_records_to_dataframe(
        records.iter().map(|r| (*r).clone())
    );

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-103-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_5() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-104-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_6() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-105-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_7() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-106-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_8() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-107-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_9() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000269305.9");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-108-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_10() {
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
        1,
        0.95f32,
        25,
        25,
        1,
        1_000,
        ""
    );

    let structure_records: Vec<TranscriptModelStructureRecord> =
        build_transcript_model_structure_records(&transcript_model_set).collect();
    let expected_rt_set: HashSet<&str> = HashSet::from([
        "ENST00000263087.9",
        "ENST00000570791.5",
        "ENST00000333813.4",
    ]);
    let records: Vec<&TranscriptModelStructureRecord> = structure_records
        .iter()
        .filter(|r| {
            let rt_set: HashSet<&str> = r.reference_transcript_id.split(',').collect();
            rt_set == expected_rt_set
        })
        .collect();
    let df_transcript_structures_filtered: DataFrame = transcript_model_structure_records_to_dataframe(
        records.iter().map(|r| (*r).clone())
    );

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-109-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_11() {
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
        1,
        0.95f32,
        25,
        25,
        1,
        1_000,
        ""
    );

    let structure_records: Vec<TranscriptModelStructureRecord> =
        build_transcript_model_structure_records(&transcript_model_set).collect();
    let expected_rt_set: HashSet<&str> = HashSet::from([
        "ENST00000263092.11",
        "ENST00000250113.12",
        "ENST00000355530.7",
    ]);
    let records: Vec<&TranscriptModelStructureRecord> = structure_records
        .iter()
        .filter(|r| {
            let rt_set: HashSet<&str> = r.reference_transcript_id.split(',').collect();
            rt_set == expected_rt_set
        })
        .collect();
    let df_transcript_structures_filtered: DataFrame = transcript_model_structure_records_to_dataframe(
        records.iter().map(|r| (*r).clone())
    );

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-110-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();

        assert!(df_results.unwrap().height() >= 1);
    }
}

#[test]
fn test_transcript_model_set_12() {
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
        1,
        1_000,
        ""
    );

    // Compare against the ground truth
    let mut df_transcript_structures: DataFrame = transcript_model_structure_records_to_dataframe(
        build_transcript_model_structure_records(&transcript_model_set)
    );
    let mask = df_transcript_structures
        .column("reference_transcript_id")
        .unwrap()
        .str()
        .unwrap()
        .equal("ENST00000254719.10");
    let df_transcript_structures_filtered = df_transcript_structures.filter(&mask).unwrap();

    let tsv_path = Path::new("src/tests/data/tsv/ground_truth/rna-111-tumor_transcript_structure_ground_truth.tsv");
    let tsv_full_path = fs::canonicalize(tsv_path).unwrap();
    let file = File::open(tsv_full_path);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());
    for result in reader.records() {
        let record = result.unwrap();
        let index: u32 = record[2].parse::<u32>().expect("Failed to convert &str to usize");
        let start: u32 = record[3].parse::<u32>().expect("Failed to convert &str to usize");
        let end: u32 = record[4].parse::<u32>().expect("Failed to convert &str to usize");
        let sequence: String = record[5].to_string().to_uppercase();
        let record_type: String = record[6].to_string();
        let kind: String = record[7].to_string();
        let context: String = record[8].to_string();
        let chromosome_1: String = record[9].to_string();
        let position_1: u32 = record[10].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_1: String = record[11].to_string();
        let strand_1: String = record[12].to_string();
        let chromosome_2: String = record[13].to_string();
        let position_2: u32 = record[14].parse::<u32>().expect("Failed to convert &str to usize");
        let operation_2: String = record[15].to_string();
        let strand_2: String = record[16].to_string();
        let gene_id_1: String = record[17].to_string();
        let transcript_id_1: String = record[18].to_string();
        let exon_id_1: String = record[19].to_string();
        let gene_id_2: String = record[20].to_string();
        let transcript_id_2: String = record[21].to_string();
        let exon_id_2: String = record[22].to_string();
        let skipped: String = record[23].to_string();

        let df_results: PolarsResult<DataFrame> = df_transcript_structures_filtered
            .clone()
            .lazy()
            .filter(col("read_start").eq(lit(start)))
            .filter(col("read_end").eq(lit(end)))
            .filter(col("sequence").eq(lit(sequence)))
            .filter(col("type").eq(lit(record_type)))
            .filter(col("kind").eq(lit(kind)))
            .filter(col("context").eq(lit(context)))
            .filter(col("chromosome_1").eq(lit(chromosome_1)))
            .filter(col("position_1").eq(lit(position_1)))
            .filter(col("operation_1").eq(lit(operation_1)))
            .filter(col("strand_1").eq(lit(strand_1)))
            .filter(col("chromosome_2").eq(lit(chromosome_2)))
            .filter(col("position_2").eq(lit(position_2)))
            .filter(col("operation_2").eq(lit(operation_2)))
            .filter(col("strand_2").eq(lit(strand_2)))
            .filter(col("gene_id_1").eq(lit(gene_id_1)))
            .filter(col("transcript_id_1").eq(lit(transcript_id_1)))
            .filter(col("exon_id_1").eq(lit(exon_id_1)))
            .filter(col("gene_id_2").eq(lit(gene_id_2)))
            .filter(col("transcript_id_2").eq(lit(transcript_id_2)))
            .filter(col("exon_id_2").eq(lit(exon_id_2)))
            .filter(col("skipped").eq(lit(skipped)))
            .collect();
        
        assert!(df_results.unwrap().height() >= 1);
    }
}

