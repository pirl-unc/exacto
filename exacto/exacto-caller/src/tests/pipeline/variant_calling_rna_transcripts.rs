use csv::ReaderBuilder;
use exacto_core::prelude::*;
use noodles_bam as bam;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::io::builders::{build_assembled_transcript_alignment_records, build_assembled_transcript_variant_records};

use super::*;


#[test]
fn scga_mini_rna_001_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-001_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-001-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-001-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_002_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-002-tumor_minimap2_sorted.bam");
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-002-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-002_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-002-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-002-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_003_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-003-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-003_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-003-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-003-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_004_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-004-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-004_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-004-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-004-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_005_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-005-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-005_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-005-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-005-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_006_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-006-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-006_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-006-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // minimap2 places the unannotated junction two bases to the left (7712249/7714999) and writes one
    // inserted and one matched base before it (the input alignment, not call-rna-vars).
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-006-tumor-1", "1884", "1885"),
        ("scga-mini-rna-006-tumor-1", "1885", "2136"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-006-tumor-1", "1885", "1885", "insertion", "7714999", "7715000"),
        ("scga-mini-rna-006-tumor-1", "1886", "1886", "match", "7714999", "7714999"),
        ("scga-mini-rna-006-tumor-1", "1886", "1887", "splicing", "7712249", "7714999"),
        ("scga-mini-rna-006-tumor-1", "1887", "2136", "match", "7712000", "7712249"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-006-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_007_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-007-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 3);

    // Every model alignment row against the ground truth written by create_rna-007_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-007-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // ASPA-WSCD1: the three bases AGG at the junction fit the genome on both sides, so the breakpoint
    // spans them (read 794-798) where the truth joins 796 to 797.
    // Wild-type WSCD1: call-rna-vars also matches ENST00000573619.1 and lists four of its exons
    // (chr17:5772926-5990021) as skipped by the splice 6070652/6080371, outside the junction.
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-007-tumor-1", "689", "796"),
        ("scga-mini-rna-007-tumor-1", "796", "797"),
        ("scga-mini-rna-007-tumor-1", "797", "911"),
        ("scga-mini-rna-007-tumor-2", "273", "274"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-007-tumor-1", "689", "794", "match", "3489235", "3489340"),
        ("scga-mini-rna-007-tumor-1", "794", "798", "breakpoint", "3489340", "6087991"),
        ("scga-mini-rna-007-tumor-1", "798", "911", "match", "6087991", "6088104"),
        ("scga-mini-rna-007-tumor-2", "273", "274", "splicing", "6070652", "6080371"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, with the positions within 2 bases: the shared bases AGG at the junction move each breakpoint by up to two bases.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-007-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && (call.position_1 as i64 - record[2].parse::<i64>().unwrap()).abs() <= 2
                    && (call.position_2 as i64 - record[6].parse::<i64>().unwrap()).abs() <= 2
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_008_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-008-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-008_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-008-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-008-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_009_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-009-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-009_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-009-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-009-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_010_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-010-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-010_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-010-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-010-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_011_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-011-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-011_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-011-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // The bases CTG (read 318-320) and ACCT (426-429) fit the genome on both sides of the two joins
    // around exon 10, so each breakpoint spans them; the first, a canonical splice in the transcript,
    // is typed backsplicing because every join of a multi-lap read reads as a back-splice.
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-011-tumor-1", "247", "320"),
        ("scga-mini-rna-011-tumor-1", "320", "321"),
        ("scga-mini-rna-011-tumor-1", "321", "427"),
        ("scga-mini-rna-011-tumor-1", "427", "428"),
        ("scga-mini-rna-011-tumor-1", "428", "537"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-011-tumor-1", "247", "317", "match", "7673538", "7673608"),
        ("scga-mini-rna-011-tumor-1", "317", "321", "breakpoint", "7670715", "7673538"),
        ("scga-mini-rna-011-tumor-1", "321", "425", "match", "7670611", "7670715"),
        ("scga-mini-rna-011-tumor-1", "425", "430", "breakpoint", "7670611", "7674288"),
        ("scga-mini-rna-011-tumor-1", "430", "537", "match", "7674181", "7674288"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, with the positions within 2 bases: the shared bases ACCT at the back-splice move each breakpoint by two bases.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-011-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && (call.position_1 as i64 - record[2].parse::<i64>().unwrap()).abs() <= 2
                    && (call.position_2 as i64 - record[6].parse::<i64>().unwrap()).abs() <= 2
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_012_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-012-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-012_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-012-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // minimap2 writes the second exon copy, exons 9 and 10 and 40 bases of exon 11 (read 1061-1825)
    // as one 765-base insertion after a splice from exon 8; call-rna-vars writes those bases twice,
    // as the splice event's sequence (1060-1826) and as an insertion row (1061-1825).
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-012-tumor-1", "1060", "1061"),
        ("scga-mini-rna-012-tumor-1", "1061", "1244"),
        ("scga-mini-rna-012-tumor-1", "1244", "1245"),
        ("scga-mini-rna-012-tumor-1", "1245", "1357"),
        ("scga-mini-rna-012-tumor-1", "1357", "1358"),
        ("scga-mini-rna-012-tumor-1", "1358", "1467"),
        ("scga-mini-rna-012-tumor-1", "1467", "1468"),
        ("scga-mini-rna-012-tumor-1", "1468", "1604"),
        ("scga-mini-rna-012-tumor-1", "1604", "1605"),
        ("scga-mini-rna-012-tumor-1", "1605", "1678"),
        ("scga-mini-rna-012-tumor-1", "1678", "1679"),
        ("scga-mini-rna-012-tumor-1", "1679", "1785"),
        ("scga-mini-rna-012-tumor-1", "1785", "1786"),
        ("scga-mini-rna-012-tumor-1", "1786", "3055"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-012-tumor-1", "1060", "1826", "splicing", "7669650", "7673701"),
        ("scga-mini-rna-012-tumor-1", "1061", "1825", "insertion", "7669650", "7669651"),
        ("scga-mini-rna-012-tumor-1", "1826", "3055", "match", "7668421", "7669650"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    // minimap2 writes the duplicated exons as an insertion (above), so no duplication is called.
    let variant_types_not_called: Vec<&str> = vec!["DUP"];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-012-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_013_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-013_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-013-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // minimap2 writes the inverted exon copy, exons 9 and 10 and 40 bases of exon 11 (read 1061-1825)
    // as one 765-base insertion after a splice from exon 8; call-rna-vars writes those bases twice,
    // as the splice event's sequence (1060-1826) and as an insertion row (1061-1825).
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-013-tumor-1", "1060", "1061"),
        ("scga-mini-rna-013-tumor-1", "1061", "1197"),
        ("scga-mini-rna-013-tumor-1", "1197", "1198"),
        ("scga-mini-rna-013-tumor-1", "1198", "1307"),
        ("scga-mini-rna-013-tumor-1", "1307", "1308"),
        ("scga-mini-rna-013-tumor-1", "1308", "1420"),
        ("scga-mini-rna-013-tumor-1", "1420", "1421"),
        ("scga-mini-rna-013-tumor-1", "1421", "1604"),
        ("scga-mini-rna-013-tumor-1", "1604", "1605"),
        ("scga-mini-rna-013-tumor-1", "1605", "1678"),
        ("scga-mini-rna-013-tumor-1", "1678", "1679"),
        ("scga-mini-rna-013-tumor-1", "1679", "1785"),
        ("scga-mini-rna-013-tumor-1", "1785", "1786"),
        ("scga-mini-rna-013-tumor-1", "1786", "3055"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-013-tumor-1", "1060", "1826", "splicing", "7669650", "7673701"),
        ("scga-mini-rna-013-tumor-1", "1061", "1825", "insertion", "7669650", "7669651"),
        ("scga-mini-rna-013-tumor-1", "1826", "3055", "match", "7668421", "7669650"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    // minimap2 writes the inverted exon copy as an insertion (above), so neither breakend is called.
    let variant_types_not_called: Vec<&str> = vec!["BND"];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-013-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_014_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-014-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 2);

    // Every model alignment row against the ground truth written by create_rna-014_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-014-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, at the same positions.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-014-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && call.position_1 == record[2].parse::<u32>().unwrap()
                    && call.position_2 == record[6].parse::<u32>().unwrap()
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_015_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 4);

    // Every model alignment row against the ground truth written by create_rna-015_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-015-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // ASPA-WSCD1: the three bases AGG at the junction fit the genome on both sides, so the breakpoint
    // spans them (read 794-798) where the truth joins 796 to 797; WSCD1-ACAP1 likewise with AG (1541-1544).
    // Wild-type WSCD1: call-rna-vars also matches ENST00000573619.1 and lists four of its exons
    // (chr17:5772926-5990021) as skipped by the splice 6070652/6080371, outside the junction.
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-015-tumor-1", "689", "796"),
        ("scga-mini-rna-015-tumor-1", "796", "797"),
        ("scga-mini-rna-015-tumor-1", "797", "911"),
        ("scga-mini-rna-015-tumor-1", "1379", "1543"),
        ("scga-mini-rna-015-tumor-1", "1543", "1544"),
        ("scga-mini-rna-015-tumor-2", "273", "274"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-015-tumor-1", "689", "794", "match", "3489235", "3489340"),
        ("scga-mini-rna-015-tumor-1", "794", "798", "breakpoint", "3489340", "6087991"),
        ("scga-mini-rna-015-tumor-1", "798", "911", "match", "6087991", "6088104"),
        ("scga-mini-rna-015-tumor-1", "1379", "1541", "match", "6110771", "6110933"),
        ("scga-mini-rna-015-tumor-1", "1541", "1544", "breakpoint", "6110933", "7341948"),
        ("scga-mini-rna-015-tumor-2", "273", "274", "splicing", "6070652", "6080371"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, with the positions within 2 bases: the shared bases AGG and AG at the junctions move each breakpoint by up to two bases.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-015-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && (call.position_1 as i64 - record[2].parse::<i64>().unwrap()).abs() <= 2
                    && (call.position_2 as i64 - record[6].parse::<i64>().unwrap()).abs() <= 2
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}

#[test]
fn scga_mini_rna_016_identify_rna_transcript_variants_returns_matches() {
    let data_dir = Path::new(env!("EXACTO_TEST_DATA"));
    let bam_file = data_dir.join("simulation/ground_truth/scga-mini-rna-016-tumor_minimap2_sorted.bam");
    let reference_genome_fasta_file = data_dir.join("references/hg38_chr17-18.fa.gz");
    let gencode_gtf_file = data_dir.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gene_annotator = Gencode::new(
        gencode_gtf_file.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file.to_str().unwrap(),
        reference_genome_fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        1,
        ""
    );

    // The input is the simulated transcripts themselves, aligned: one model per transcript.
    assert_eq!(transcript_model_set.transcript_models.len(), 3);

    // Every model alignment row against the ground truth written by create_rna-016_fasta_file.py.
    // Rows compare gene and exon ids, not transcript ids: an exon shared by several isoforms may be
    // named under any of them.
    let rows: Vec<Vec<String>> = build_assembled_transcript_alignment_records(&transcript_model_set)
        .map(|record| vec![
            record.assembled_transcript_name.to_string(), record.read_start.to_string(), record.read_end.to_string(),
            record.sequence.to_string(), record.record_type.to_string(), record.kind.to_string(), record.context.to_string(),
            record.chromosome_1.to_string(), record.position_1.to_string(), record.operation_1.to_string(), record.strand_1.to_string(),
            record.chromosome_2.to_string(), record.position_2.to_string(), record.operation_2.to_string(), record.strand_2.to_string(),
            record.reference_gene_id_1.to_string(), record.reference_exon_id_1.to_string(),
            record.reference_gene_id_2.to_string(), record.reference_exon_id_2.to_string(), record.skipped.to_string()
        ])
        .collect();
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-016-tumor_transcript_model_alignments_ground_truth.tsv"))
        .unwrap();
    let truth_rows: Vec<Vec<String>> = reader
        .records()
        .map(|result| {
            let record = result.unwrap();
            [0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 23, 24].iter().map(|&i| record[i].to_string()).collect()
        })
        .collect();

    // Known differences: truth rows written otherwise, by (transcript, read_start, read_end), and the
    // rows written instead, by (transcript, read_start, read_end, kind, position_1, position_2).
    // The last inserted base G equals the intron base before WSCD1 exon 3 (6087989), so the breakpoint
    // ends one base early and the G is placed in the intron.
    // Wild-type WSCD1: call-rna-vars also matches ENST00000573619.1 and lists four of its exons
    // (chr17:5772926-5990021) as skipped by the splice 6070652/6080371, outside the junction.
    let truth_rows_written_otherwise: Vec<(&str, &str, &str)> = vec![
        ("scga-mini-rna-016-tumor-1", "796", "809"),
        ("scga-mini-rna-016-tumor-2", "273", "274"),
    ];
    let rows_written_instead: Vec<(&str, &str, &str, &str, &str, &str)> = vec![
        ("scga-mini-rna-016-tumor-1", "796", "808", "breakpoint", "3489342", "6087989"),
        ("scga-mini-rna-016-tumor-1", "808", "808", "match", "6087989", "6087989"),
        ("scga-mini-rna-016-tumor-2", "273", "274", "splicing", "6070652", "6080371"),
    ];
    for truth_row in truth_rows.iter() {
        let expected: usize = if truth_rows_written_otherwise.contains(&(truth_row[0].as_str(), truth_row[1].as_str(), truth_row[2].as_str())) { 0 } else { 1 };
        assert_eq!(rows.iter().filter(|row| *row == truth_row).count(), expected, "ground truth row {:?}", truth_row);
    }
    for row in rows.iter().filter(|row| !truth_rows.contains(row)) {
        assert!(
            rows_written_instead.contains(&(row[0].as_str(), row[1].as_str(), row[2].as_str(), row[5].as_str(), row[8].as_str(), row[12].as_str())),
            "row not in the ground truth {:?}", row
        );
    }
    assert_eq!(rows.len(), truth_rows.len() - truth_rows_written_otherwise.len() + rows_written_instead.len());

    // Every variant ground truth row is exactly one call: same chromosomes, strands, operations and
    // type, with the positions within 1 base: the last inserted base matches the intron, which moves the second breakpoint by one base.
    let calls: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&transcript_model_set, None, 20).collect();
    let variant_types_not_called: Vec<&str> = vec![];
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(data_dir.join("simulation/ground_truth/scga-mini-rna-016-tumor_ground_truth.tsv"))
        .unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let num_matches: usize = calls
            .iter()
            .filter(|call| {
                *call.chromosome_1 == record[1]
                    && *call.strand_1 == record[3]
                    && *call.operation_1 == record[4]
                    && *call.chromosome_2 == record[5]
                    && *call.strand_2 == record[7]
                    && *call.operation_2 == record[8]
                    && *call.variant_type == record[10]
                    && (call.position_1 as i64 - record[2].parse::<i64>().unwrap()).abs() <= 1
                    && (call.position_2 as i64 - record[6].parse::<i64>().unwrap()).abs() <= 1
            })
            .count();
        let expected: usize = if variant_types_not_called.contains(&&record[10]) { 0 } else { 1 };
        assert_eq!(num_matches, expected, "variant ground truth row {:?}", record);
    }
}


/// Every record of scga-mini-rna-015 (its four simulated transcripts, aligned) given a cs tag with an unknown reference base, so every read fails.
/// A read that fails is left out, but when most reads fail the input is at fault: the run
/// panics with a count of them instead of returning a set that looks complete.
#[test]
#[should_panic(expected = "4 of the 4 reads of")]
fn scga_mini_rna_015_identify_rna_transcript_variants_panics_when_most_reads_fail() {
    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use tempfile::TempDir;

    let bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-015-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let reference_genome_fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator = Gencode::new(
        gencode_gtf_full_path.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let temp_dir = TempDir::new().unwrap();
    let broken_bam_file: String = temp_dir.path().join("broken_cs.bam").to_str().unwrap().to_string();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&broken_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        record.data_mut().insert(Tag::from([b'c', b's']), Value::String("*zz".into()));
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);

    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    identify_rna_transcript_variants(
        &broken_bam_file,
        reference_genome_fasta_file,
        &gene_annotator,
        &options,
        1,
        temp_dir.path().to_str().unwrap()
    );
}


/// Two runs of one input give the same tables, row for row: the models are held in read order.
#[test]
fn scga_mini_rna_001_identify_rna_transcript_variants_gives_the_same_tables_on_two_runs() {
    let bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let reference_genome_fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator = Gencode::new(
        gencode_gtf_full_path.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();

    let mut tables: Vec<(Vec<AssembledTranscriptRecord>, Vec<AssembledTranscriptVariantRecord>, Vec<AssembledTranscriptFilterStatusRecord>)> = Vec::new();
    for _ in 0..2 {
        let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
            bam_file,
            reference_genome_fasta_file,
            &gene_annotator,
            &options,
            2,
            ""
        );
        tables.push((
            build_assembled_transcript_records(&transcript_model_set).collect(),
            build_assembled_transcript_variant_records(&transcript_model_set, None, 0).collect(),
            build_assembled_transcript_filter_status_records(&transcript_model_set).collect()
        ));
    }

    assert!(tables[0].0.len() > 100);
    assert!(!tables[0].1.is_empty());
    assert_eq!(tables[1], tables[0]);
}

/// The simulated transcripts of scga-mini-rna-015 (two of them fusions split over two records) with
/// their supplementary records hard-clipped, as minimap2 writes them without -Y, give the same
/// transcript models as with soft clips. The split reads were the ones dropped before.
#[test]
fn scga_mini_rna_015_identify_rna_transcript_variants_reads_hard_clipped_supplementary_records() {
    use noodles_sam::alignment::io::Write;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use tempfile::TempDir;

    let bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-015-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let reference_genome_fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator = Gencode::new(
        gencode_gtf_full_path.to_str().unwrap(),
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let temp_dir = TempDir::new().unwrap();
    let hard_clipped_bam_file: String = temp_dir.path().join("hard_clipped.bam").to_str().unwrap().to_string();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&hard_clipped_bam_file).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        if record.flags().is_supplementary() {
            let operations: Vec<Op> = record.cigar().as_ref().to_vec();
            let leading: usize = if operations[0].kind() == Kind::SoftClip { operations[0].len() } else { 0 };
            let trailing: usize = if operations[operations.len() - 1].kind() == Kind::SoftClip { operations[operations.len() - 1].len() } else { 0 };
            let length: usize = record.sequence().len();
            let sequence: Vec<u8> = record.sequence().as_ref()[leading..length - trailing].to_vec();
            // The simulated transcripts were aligned from a FASTA: no base qualities to clip.
            let quality_scores: Vec<u8> = record.quality_scores().as_ref().get(leading..length - trailing).unwrap_or_default().to_vec();
            *record.sequence_mut() = sequence.into();
            *record.quality_scores_mut() = quality_scores.into();
            *record.cigar_mut() = operations
                .into_iter()
                .map(|op| if op.kind() == Kind::SoftClip { Op::new(Kind::HardClip, op.len()) } else { op })
                .collect();
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);

    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let soft_clipped: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file, reference_genome_fasta_file, &gene_annotator, &options, 1, temp_dir.path().to_str().unwrap()
    );
    let hard_clipped: TranscriptModelSet = identify_rna_transcript_variants(
        &hard_clipped_bam_file, reference_genome_fasta_file, &gene_annotator, &options, 1, temp_dir.path().to_str().unwrap()
    );

    let soft_clipped_rows: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&soft_clipped, None, 0).collect();
    let hard_clipped_rows: Vec<AssembledTranscriptVariantRecord> = build_assembled_transcript_variant_records(&hard_clipped, None, 0).collect();
    assert_eq!(soft_clipped.get_size(), 4);
    assert_eq!(hard_clipped.get_size(), 4);
    assert_eq!(hard_clipped_rows, soft_clipped_rows);
}

/// The nonsense-mediated decay prediction of every open reading frame is written as a table.
#[test]
fn scga_mini_rna_001_identify_rna_transcript_variants_writes_the_nmd_predictions() {
    let bam_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam")).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let reference_genome_fasta_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz")).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_full_path = fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz")).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gencode_gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();

    let transcript_model_set: TranscriptModelSet = identify_rna_transcript_variants(
        bam_file,
        reference_genome_fasta_file,
        &gene_annotator,
        &options,
        2,
        ""
    );

    let rows: Vec<(Box<str>, u32, u32, bool)> = build_assembled_transcript_nonsense_mediated_decay_records(&transcript_model_set)
        .map(|record| (record.assembled_transcript_name, record.orf_start, record.orf_end, record.nmd_predicted))
        .collect();
    let start_codons: HashSet<&str> = options.annotation.start_codons.iter().map(|codon| codon.as_str()).collect();
    let nmd_predictor: NonsenseMediatedDecayPredictor = NonsenseMediatedDecayPredictor {
        translation_strategy: &options.annotation.translation_strategy,
        start_codons: &start_codons,
        distance_threshold: options.annotation.nmd_distance_threshold
    };
    let mut expected: Vec<(Box<str>, u32, u32, bool)> = Vec::new();
    for transcript_model in transcript_model_set.transcript_models.iter() {
        let read_name: Box<str> = transcript_model_set.read_names_map.get_by_right(&transcript_model.get_read_id()).unwrap().clone();
        for call in nmd_predictor.predict(transcript_model).iter() {
            let predicted: bool = matches!(call.verdict, NonsenseMediatedDecayVerdict::Predicted { .. });
            expected.push((read_name.clone(), call.orf_start, call.orf_end, predicted));
        }
    }

    assert_eq!(transcript_model_set.transcript_models.len(), 2);
    assert!(!rows.is_empty());
    assert_eq!(rows, expected);
}
