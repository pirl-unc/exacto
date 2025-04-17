use std::fs;
use std::fs::File;
use std::path::Path;
use bimap::BiMap;
use csv::ReaderBuilder;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::{col, lit};
use exacto_util::prelude::{Gencode, GeneAnnotator};
use crate::prelude::{create_chromosome_names_map, create_read_names_map, identify_variant_transcripts, TranscriptModelSet};


// #[test]
// fn test_playground() {
//     let bam_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/test2.bam";
//     let bam_bai_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/experiments/simulate_transcriptome/bam/test2.bam.bai";
//     let reference_genome_fasta_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/mm39.fa";
//     let gencode_gtf_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/gencode.vM36.annotation.gtf.gz";
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
//         25,
//         25f32,
//         2
//     );
// }