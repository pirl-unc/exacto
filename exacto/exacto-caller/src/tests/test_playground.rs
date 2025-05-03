use std::fs;
use std::fs::File;
use std::path::Path;
use bimap::BiMap;
use csv::ReaderBuilder;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::{col, lit};
use exacto_util::prelude::{Gencode, GeneAnnotator};
use crate::prelude::{create_chromosome_names_map, create_read_names_map, identify_variant_transcripts, ReferenceTranscriptScoringMethod, TranscriptModelSet};

// use crate::algorithms::variant_calling_rna::svd_l2_match;
// use nalgebra::{DMatrix, DVector, SVD};


// #[test]
// fn test_playground() {
//     let R = DMatrix::from_row_slice(3, 11, &vec![
//         1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
//         1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
//         1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0
//     ]);
// 
//     // Step 2: Create vector B (length 11)
//     let B = DVector::from_vec(vec![
//         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0
//     ]);
// 
//     // Step 3: Compute results
//     let (best_idx, distances) = svd_l2_match(&R, &B, 2).unwrap();
// 
// }

// #[test]
// fn test_playground() {
//     // let bam_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/mm39_protein_coding_transcripts_0percent_degradation_minimap2_mdtagged_sorted.bam";
//     // let bam_bai_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/mm39_protein_coding_transcripts_0percent_degradation_minimap2_mdtagged_sorted.bam.bai";
//     // let reference_genome_fasta_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/mm39.fa";
//     // let gencode_gtf_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/gencode.vM36.annotation.gtf.gz";
// 
//     let bam_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/pmch.bam";
//     let bam_bai_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/pmch.bam.bai";
//     let reference_genome_fasta_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/mm39.fa";
//     let gencode_gtf_file: &str = "/Users/leework/Desktop/temp/20250423/pmch.gtf";
// 
// 
//     let gene_annotator = Gencode::new(
//         gencode_gtf_file,
//         "mm39"
//     );
// 
//     let mut transcript_model_set: TranscriptModelSet = identify_variant_transcripts(
//         bam_file,
//         bam_bai_file,
//         reference_genome_fasta_file,
//         &gene_annotator,
//         ReferenceTranscriptScoringMethods::SVD,
//         25,
//         25f32,
//         1
//     );
// }