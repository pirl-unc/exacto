// use noodles_bam as bam;
// use std::collections::HashMap;
// use std::fs;
// use std::path::Path;
// use exacto_util::prelude::*;
// use crate::common::bam::*;
// use crate::structs::alignment::Alignment;
// use crate::structs::variant_record::VariantRecord;
// use crate::algorithms::variant_calling_dna::{identify_dna_variants,identify_case_specific_dna_variants};
// use crate::algorithms::variant_calling_rna::identify_variant_transcripts;
// use crate::structs::variant_call::VariantCall;
// use crate::structs::variant_call_set::VariantCallSet;


// #[test]
// fn test_playground() {
    // let bam_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_giab/alignment/dna/pacbio/hg001_dna_pacbio-nist_minimap2_mdtagged_sorted_chr21_chr22.bam");
    // let bam_full_path = fs::canonicalize(bam_path).unwrap();
    // let bam_file: &str = bam_full_path.to_str().unwrap();
    // let bam_bai_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_giab/alignment/dna/pacbio/hg001_dna_pacbio-nist_minimap2_mdtagged_sorted_chr21_chr22.bam.bai");
    // let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    // let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    // let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    // let end: usize = *chromosome_lengths.get("chr21").unwrap();
    // let read_ids_map: HashMap<Box<str>,usize> = get_read_ids_map(bam_file);
    // let chromosomes: Vec<&str> = vec!["chr11","chr13"];
    // let mut regions: HashMap<Box<str>,Vec<(usize,usize)>> = generate_regions(
    //     bam_file,
    //     &chromosomes,
    //     10_000_000
    // );
    // for (_, vec) in regions.iter_mut() {
    //     vec.sort_by(|a, b| a.0.cmp(&b.0));
    // }
    // println!("{:?}", regions);

//     let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
//         bam_file,
//         bam_bai_file,
//         "chr21",
//         1,
//         46709983,
//         &read_ids_map
//     );
// }

// #[test]
// fn test_playground() {
//     let bam_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/bbn963tumor_rna_pacbio-2024duke_flnc_minimap2_mdtagged_sorted_fancl_bnd.bam");
//     let bam_full_path = fs::canonicalize(bam_path).unwrap();
//     let bam_file: &str = bam_full_path.to_str().unwrap();
//     let bam_bai_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/exacto/exacto/exacto-caller/src/tests/data/bam/bbn963tumor_rna_pacbio-2024duke_flnc_minimap2_mdtagged_sorted_fancl_bnd.bam.bai");
//     let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
//     let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
//     let records_map: HashMap<Box<str>,Vec<bam::Record>> = fetch_bam_records(
//         bam_file,
//         bam_bai_file,
//         1,
//         &vec!["chr11"]
//     );
//     let chromosomes_ids_names: HashMap<usize,Box<str>> = get_chromosome_ids_names(bam_file);
//     let read_id: &str = "m84165_240529_114357_s1/153555735/ccs/12853_14228";
//     let read_sequence: Box<str> = get_original_read_sequence(records_map.get(read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
//     let quality_scores: Vec<u8> = get_original_base_quality_scores(records_map.get(read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
//     let alignment: Alignment = Alignment::new(
//         read_id,
//         &*read_sequence,
//         &quality_scores,
//         records_map.get(read_id).unwrap()
//     );
//     for alignment_record in alignment.alignment_records.iter() {
//         println!("left_softclipping.0: {}", get_left_softclipping(&alignment_record.record).0);
//         println!("left_softclipping.1: {}", get_left_softclipping(&alignment_record.record).1);
//         println!("right_softclipping.0: {}", get_right_softclipping(&alignment_record.record).0);
//         println!("right_softclipping.1: {}", get_right_softclipping(&alignment_record.record).1);
//         println!("alignment start: {}", get_alignment_start_position(&alignment_record.record));
//         println!("alignment_record.read_start: {}", alignment_record.read_start);
//         println!("alignment_record.read_end: {}", alignment_record.read_end);
//         println!("alignment sequence: {}", get_aligned_sequence_from_cigar(&alignment_record.record));
//     }

    // let variant_records: Vec<VariantRecord> = alignment.identify_variant_records_in_softclipping(
    //     30,
    //     &chromosomes_ids_names
    // );
    // println!("{} variant records", variant_records.len());
    // for variant_record in variant_records.iter() {
    //     println!("{}", variant_record.get_sequence_operation_string());
    // }
    // println!("variant_records: {:?}", variant_records);
    // assert!(alignment.get_alignment_records_count() == 2);
// }




// #[test]
// fn test_playground() {
//     let case_bam_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/dna/pacbio/minimap2/mm39_chr1-19_chrX/bbn963tumor_dna_pacbio-2024duke_minimap2_mdtagged_sorted_fancl.bam");
//     let case_bam_full_path = fs::canonicalize(case_bam_path).unwrap();
//     let case_bam_file: &str = case_bam_full_path.to_str().unwrap();
//     let case_bam_bai_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/dna/pacbio/minimap2/mm39_chr1-19_chrX/bbn963tumor_dna_pacbio-2024duke_minimap2_mdtagged_sorted_fancl.bam.bai");
//     let case_bam_bai_full_path = fs::canonicalize(case_bam_bai_path).unwrap();
//     let case_bam_bai_file: &str = case_bam_bai_full_path.to_str().unwrap();
//
//     // let case_bam_path = Path::new("/Users/leework/Desktop/sandbox/20250204/test_sorted.bam");
//     // let case_bam_full_path = fs::canonicalize(case_bam_path).unwrap();
//     // let case_bam_file: &str = case_bam_full_path.to_str().unwrap();
//     // let case_bam_bai_path = Path::new("/Users/leework/Desktop/sandbox/20250204/test_sorted.bam.bai");
//     // let case_bam_bai_full_path = fs::canonicalize(case_bam_bai_path).unwrap();
//     // let case_bam_bai_file: &str = case_bam_bai_full_path.to_str().unwrap();
//
//     let control_bam_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/dna/pacbio/minimap2/mm39_chr1-19_chrX/bl6_dna_pacbio-2024duke_minimap2_mdtagged_sorted_fancl_2.bam");
//     let control_bam_full_path = fs::canonicalize(control_bam_path).unwrap();
//     let control_bam_file: &str = control_bam_full_path.to_str().unwrap();
//     let control_bam_bai_path = Path::new("/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/dna/pacbio/minimap2/mm39_chr1-19_chrX/bl6_dna_pacbio-2024duke_minimap2_mdtagged_sorted_fancl_2.bam.bai");
//     let control_bam_bai_full_path = fs::canonicalize(control_bam_bai_path).unwrap();
//     let control_bam_bai_file: &str = control_bam_bai_full_path.to_str().unwrap();
//
//     let control_bam_files: Vec<&str> = vec![control_bam_file];
//     let control_bam_bai_files: Vec<&str> = vec![control_bam_bai_file];
//
//     let chromosomes: Vec<&str> = vec!["chr11"];
//
//     // let variant_call_set: VariantCallSet = identify_dna_variants(
//     //     case_bam_file,
//     //     case_bam_bai_file,
//     //     1,
//     //     20,
//     //     30f32,
//     //     0.5f32,
//     //     0.5,
//     //     2000,
//     //     1000,
//     //     1000,
//     //     4,
//     //     chromosomes,
//     //     ""
//     // );
//
//     let variant_call_set: VariantCallSet = identify_case_specific_dna_variants(
//         case_bam_file,
//         case_bam_bai_file,
//         control_bam_files,
//         control_bam_bai_files,
//         3,
//         25,
//         20f32,
//         0.5,
//         0.5,
//         2000,
//         1000,
//         1000,
//         true,
//         2,
//         chromosomes,
//         ""
//     );
//
//     for variant_call in variant_call_set.variant_calls.iter() {
//         if variant_call.get_consensus_record().0.get_variant_type() == VariantTypes::BREAKPOINT {
//             println!("{:?}", variant_call.get_consensus_record().0);
//         }
//     }
// }


// #[test]
// fn test_playground() {
//     let bam_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/rna/pacbio/minimap2/mm39_chr1-19_chrX/bbn963tumor_rna_pacbio-2024duke_flnc_rnabloom2_filtered_minimap2_mdtagged_sorted.bam";
//     let bam_bai_file: &str = "/Users/leework/Documents/Research/projects/project_exacto/data/processed/samples_mouse/alignment/rna/pacbio/minimap2/mm39_chr1-19_chrX/bbn963tumor_rna_pacbio-2024duke_flnc_rnabloom2_filtered_minimap2_mdtagged_sorted.bam.bai";
//     let refernce_fasta_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/mm39.fa";
//     let gencode_file: &str = "/Users/leework/Documents/Research/projects/seqdata/references/gencode.vM36.gene_annotation.gtf";
//     let gencode: Gencode = Gencode::new(gencode_file,"mm39", 2);
//     identify_variant_transcripts(
//         bam_file,
//         bam_bai_file,
//         refernce_fasta_file,
//         &gencode,
//         25,
//         25f32,
//         2
//     );
// }