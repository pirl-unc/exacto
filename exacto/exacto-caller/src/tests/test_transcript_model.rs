use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use crate::prelude::*;


#[test]
fn test_transcript_model_1() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_507476_774164/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    // ENST00000610538.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000610538.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 4);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668404);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673207);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673266);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(2).unwrap().get_standardized_sequence(), "A");
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7687482);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000620739.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000620739.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668402);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(1).unwrap().get_standardized_sequence(), "A");
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687491);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687538);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::ExonTruncation);

    // ENST00000455263.6
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000455263.6".into()]).unwrap();
    assert_eq!(variant_records.len(), 4);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668404);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673207);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673266);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(2).unwrap().get_standardized_sequence(), "A");
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7687482);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000619485.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000619485.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "A");
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7676620);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7676622);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::IntronRetention);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687488);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000269305.9
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "A");

    // ENST00000445888.6
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000445888.6".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7674224);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7674226);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "A");
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7676620);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7676622);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::IntronRetention);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687488);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::UTRExtension);
}

#[test]
fn test_transcript_model_2() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_822724_603243/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    // ENST00000445888.6
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000445888.6".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7676400);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7676401);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::Insertion);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "GGGGGTTTTT");

    // ENST00000610538.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000610538.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 4);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668404);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673207);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673266);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7676400);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7676401);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::Insertion);
    assert_eq!(variant_records.get(2).unwrap().get_standardized_sequence(), "GGGGGTTTTT");
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7687482);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000269305.9
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7676400);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7676401);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::Insertion);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "GGGGGTTTTT");

    // ENST00000620739.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000620739.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668402);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7676400);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7676401);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::Insertion);
    assert_eq!(variant_records.get(1).unwrap().get_standardized_sequence(), "GGGGGTTTTT");
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687491);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687538);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::ExonTruncation);
}

#[test]
fn test_transcript_model_3() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_264855_304921/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    // ENST00000455263.6
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000455263.6".into()]).unwrap();
    assert_eq!(variant_records.len(), 4);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668404);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673207);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673266);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::Deletion);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7687482);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000620739.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000620739.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668402);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::Deletion);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687491);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687538);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::ExonTruncation);

    // ENST00000445888.6
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000445888.6".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::Deletion);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7676620);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7676622);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::IntronRetention);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687488);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000269305.9
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::Deletion);

    // ENST00000610538.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000610538.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 4);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7668404);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7668420);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7673207);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7673266);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::Deletion);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7687482);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::UTRExtension);

    // ENST00000619485.4
    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000619485.4".into()]).unwrap();
    assert_eq!(variant_records.len(), 3);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7673750);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7673761);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::Deletion);

    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7676620);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7676622);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::IntronRetention);

    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7687488);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7687490);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::UTRExtension);
}

#[test]
fn test_transcript_model_4() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_535544_475898/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000698746.1".into(), "ENST00000570791.5".into()]).unwrap();
    assert_eq!(variant_records.len(), 10);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7701730);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7701732);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::SingleNucleotideVariant);

    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7701731);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7727201);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::FusionGene);

    assert_eq!(variant_records.get(2).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(2).unwrap().get_position_1(), 7701732);
    assert_eq!(variant_records.get(2).unwrap().get_position_2(), 7701789);
    assert_eq!(variant_records.get(2).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(2).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(2).unwrap().get_standardized_sequence(), "GTGCGTGTTTTTTCCACGGCCCGGCCTGGCCGAGACTGCGAGGTCCGAGCCACATTTG");

    assert_eq!(variant_records.get(3).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(3).unwrap().get_position_1(), 7702344);
    assert_eq!(variant_records.get(3).unwrap().get_position_2(), 7702552);
    assert_eq!(variant_records.get(3).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(3).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(3).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(3).unwrap().get_standardized_sequence(), "CAAAAAAGCAGGGCCAGAGCGGCATCATCTCCTGCATAGCCTTCAGCCCAGCCCAGCCCCTCTATGCCTGTGGCTCCTACGGCCGCTCCCTGGGTCTGTATGCCTGGGATGATGGCTCCCCTCTCGCCTTGCTGGGAGGGCACCAAGGGGGCATCACCCACCTCTGCTTTCATCCCGATGGCAACCGCTTCTTCTCAGGAGCCCGCAAG");

    assert_eq!(variant_records.get(4).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_position_1(), 7702743);
    assert_eq!(variant_records.get(4).unwrap().get_position_2(), 7702846);
    assert_eq!(variant_records.get(4).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(4).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(4).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(4).unwrap().get_standardized_sequence(), "GATGCTGAGCTCCTGTGCTGGGATCTCCGGCAGTCTGGTTACCCACTGTGGTCCCTGGGTCGAGAGGTGACCACCAATCAGCGCATCTACTTCGATCTGGACCC");

    assert_eq!(variant_records.get(5).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(5).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(5).unwrap().get_position_1(), 7703243);
    assert_eq!(variant_records.get(5).unwrap().get_position_2(), 7703499);
    assert_eq!(variant_records.get(5).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(5).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(5).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(5).unwrap().get_standardized_sequence(), "CCTGCACCCTAGCCTGCCTCTCCTGGCCACTGCCTCCGGTCAGCGTGTGTTTCCTGAGCCCACAGAGAGTGGGGACGAAGGAGAGGAGCTGGGCCTTCCCTTGCTCTCCACGCGCCACGTCCACCTTGAATGTCGGCTTCAGCTCTGGTGGTGTGGGGGGGCGCCAGACTCCAGCATCCCTGATGATCACCAGGGCGAGAAAGGGCAGGGAGGAACGGAGGGAGGTGTGGGTGAGCTGATATAAAAAGGTTTTTATG");

    assert_eq!(variant_records.get(6).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_position_1(), 7717761);
    assert_eq!(variant_records.get(6).unwrap().get_position_2(), 7717831);
    assert_eq!(variant_records.get(6).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(6).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(6).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(6).unwrap().get_standardized_sequence(), "AGGGCTGCGAGGGGCAACTTCTTAGAGTGGCCCATCGGTCGGTCTAGGGAGGGGAGGGTCAGCGTGGCAAG");

    assert_eq!(variant_records.get(7).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(7).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(7).unwrap().get_position_1(), 7719721);
    assert_eq!(variant_records.get(7).unwrap().get_position_2(), 7719900);
    assert_eq!(variant_records.get(7).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(7).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(7).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(7).unwrap().get_standardized_sequence(), "GTTTTGCCTGCACGATGTCCAGCAAAGCTGAGAAGAAGCAGCGATTGAGTGGCCGAGGAAGCTCCCAGGCAAGCTGGTCAGGGCGGGCCACTCGGGCTGCTGTGGCCACACAGGAGCAGGGGAATGCCCCGGCTGTCAGTGAGCCAGAGCTGCAGGCTGAGCTCCCCAAGGAGGAGCCTG");

    assert_eq!(variant_records.get(8).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(8).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(8).unwrap().get_position_1(), 7723628);
    assert_eq!(variant_records.get(8).unwrap().get_position_2(), 7723689);
    assert_eq!(variant_records.get(8).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(8).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(8).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(8).unwrap().get_standardized_sequence(), "AGCCACGGTTGGAGGGACCTCAAGCACAGAGTGAAGAATCAGTGGAGCCCGAGGCAGATGTG");

    assert_eq!(variant_records.get(9).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(9).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(9).unwrap().get_position_1(), 7727122);
    assert_eq!(variant_records.get(9).unwrap().get_position_2(), 7727200);
    assert_eq!(variant_records.get(9).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(9).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(9).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(9).unwrap().get_standardized_sequence(), "AAGCCCCTCTTCCTTTCCCGAGCTGCGCTGACAGGACTGGCGGATGCAGTGTGGACACAGGAGCATGATGCCATTCTGG");
}

#[test]
fn test_transcript_model_5() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_561742_839878/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7675994);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7676099);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "CGTGCAAGTCACAGACTTGGCTGTCCCAGAATGCAAGAAGCCCAGACGGAAACCGTAGCTGCCCTGGTAGGTTTTCTGGGAAGGGACAGAAGATGACAGGGGCCAG");
}

#[test]
fn test_transcript_model_6() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_124525_407996/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7676201);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7676272);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "TTCAGTGAACCATTGTTCAATATCGTCCGGGGACAGCATCAAATCATCCATTGCTTGGGACGGCAAGGGGGA");
}

#[test]
fn test_transcript_model_7() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_924107_174289/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7675994);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7676272);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "CGTGCAAGTCACAGACTTGGCTGTCCCAGAATGCAAGAAGCCCAGACGGAAACCGTAGCTGCCCTGGTAGGTTTTCTGGGAAGGGACAGAAGATGACAGGGGCCAGGAGGGGGCTGGTGCAGGGGCCGCCGGTGTAGGAGCTGCTGGTGCAGGGGCCACGGGGGGAGCAGCCTCTGGCATTCTGGGAGCTTCATCTGGACCTGGGTCTTCAGTGAACCATTGTTCAATATCGTCCGGGGACAGCATCAAATCATCCATTGCTTGGGACGGCAAGGGGGA");
}

#[test]
fn test_transcript_model_8() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_924958_759981/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 2);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7675601);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7675640);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::CrypticExon);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 7675993);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 7675993);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::IntronRetention);
}

#[test]
fn test_transcript_model_9() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_721712_133913/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000269305.9".into()]).unwrap();
    assert_eq!(variant_records.len(), 1);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 7675976);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 7675993);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Include);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::IntronRetention);
}

#[test]
fn test_transcript_model_10() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_288476_571946/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    assert_eq!(alignment_structure.get_base(837).is_embedded_insertion(), true);

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map
        .get(&vec![
            "ENST00000263087.9".into(),
            "ENST00000570791.5".into(),
            "ENST00000333813.4".into()
        ])
        .unwrap();

    assert_eq!(variant_records.len(), 34);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 1295600);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 3801187);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(0).unwrap().get_standardized_sequence(), "G"); // overlapping alignment
    assert_eq!(variant_records.get(29).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(29).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(29).unwrap().get_position_1(), 3761101);
    assert_eq!(variant_records.get(29).unwrap().get_position_2(), 7727201);
    assert_eq!(variant_records.get(29).unwrap().get_operation_1(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(29).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(29).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(29).unwrap().get_standardized_sequence(), "C"); // overlapping alignment
}

#[test]
fn test_transcript_model_11() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_175366_924183/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map
        .get(&vec![
            "ENST00000263092.11".into(),
            "ENST00000250113.12".into(),
            "ENST00000355530.7".into(),
        ])
        .unwrap();

    assert_eq!(variant_records.len(), 17);
    assert_eq!(variant_records.get(4).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(4).unwrap().get_position_1(), 2464208);
    assert_eq!(variant_records.get(4).unwrap().get_position_2(), 4433940);
    assert_eq!(variant_records.get(4).unwrap().get_operation_1(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(4).unwrap().get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(4).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(4).unwrap().get_standardized_sequence(), "");
    assert_eq!(variant_records.get(6).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(6).unwrap().get_position_1(), 4453100);
    assert_eq!(variant_records.get(6).unwrap().get_position_2(), 7603799);
    assert_eq!(variant_records.get(6).unwrap().get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(6).unwrap().get_operation_2(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(6).unwrap().get_variant_type(), &VariantType::FusionGene);
    assert_eq!(variant_records.get(6).unwrap().get_standardized_sequence(), "");
}

#[test]
fn test_transcript_model_12() {
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

    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: usize = *chromosome_lengths.get("chr17").unwrap();
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        1
    );
    let gene_annotator = Gencode::new(
        gencode_gtf_file,
        "hg38",
        "v41",
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2])),
        Some(HashSet::from(["protein_coding"])),
        Some(HashSet::from([1,2]))
    );
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
        bam_file,
        bam_bai_file,
        "chr17",
        1,
        end,
        &read_names_map,
        1
    );

    let read_name: &str = "m64012_324970_273886/1/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap().iter().collect::<Vec<_>>().as_slice());

    let mut alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap()
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        1,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &transcript_model.get_exons(),
        &gene_annotator,
        &chromosome_names_map,
        ReferenceTranscriptScoringMethod::CosineSimilarity,
        ReferenceTranscriptSelectionStrategy::TopK,
        3,
        0.9f32
    );

    let variant_records_map: &HashMap<Vec<Box<str>>, Vec<VariantRecord>> = transcript_model.identify_variants(
        &reference_transcript_matches,
        &gene_annotator,
        reference_genome_fasta_file,
        30,
        30
    );

    let variant_records: &Vec<VariantRecord> = variant_records_map.get(&vec!["ENST00000254719.10".into()]).unwrap();

    assert_eq!(variant_records.len(), 15);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(0).unwrap().get_position_1(), 1830005);
    assert_eq!(variant_records.get(0).unwrap().get_position_2(), 1830126);
    assert_eq!(variant_records.get(0).unwrap().get_operation_1(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_operation_2(), &GraphOperationType::Skip);
    assert_eq!(variant_records.get(0).unwrap().get_variant_type(), &VariantType::ExonTruncation);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_1(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_chromosome_2(), 0);
    assert_eq!(variant_records.get(1).unwrap().get_position_1(), 1842804);
    assert_eq!(variant_records.get(1).unwrap().get_position_2(), 1844684);
    assert_eq!(variant_records.get(1).unwrap().get_operation_1(), &GraphOperationType::Upstream);
    assert_eq!(variant_records.get(1).unwrap().get_operation_2(), &GraphOperationType::Downstream);
    assert_eq!(variant_records.get(1).unwrap().get_variant_type(), &VariantType::CircularRNA);
    assert_eq!(variant_records.get(1).unwrap().get_standardized_sequence(), "AGG");
}
