use bimap::BiMap;
use exacto_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use tempfile::TempDir;

use super::*;


/// A breakpoint whose junction carries 8 bases of flank homology (3 before it, 5 after it) and
/// no splice motif: a template switch at `junction_min_homology` 8, suppressed with its reason.
///
/// chrA (leaves at 30, '+'): positions 21-30 = "GGGGGGGCAT", 31-40 = "ACGTATTTTT".
/// chrB (enters at 51, '+'): positions 41-50 = "TTTTTTTCAT", 51-60 = "ACGTACCCCC".
#[test]
fn template_switch_filter_suppresses_a_breakpoint_with_strong_homology() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let gtf_file: String = directory.path().join("annotation.gtf").to_str().unwrap().to_string();
    // One gene, on a contig away from the junction: no exon edge at either breakend.
    std::fs::write(&gtf_file, concat!(
        "chrZ\tTEST\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_type \"protein_coding\"; level 1;\n",
        "chrZ\tTEST\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_type \"protein_coding\"; level 1;\n",
        "chrZ\tTEST\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E1\"; exon_number 1; level 1;\n"
    )).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(&gtf_file, "test", "v0");
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);
    let transcript_models_map: HashMap<ReadID, Arc<TranscriptModel>> = HashMap::new();

    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for read_id in [1, 2] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, 30, Strand::Forward, GraphOperationType::Downstream,
                1, 51, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::Breakpoint
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);
    let filter: TemplateSwitchFilter<Gencode> = TemplateSwitchFilter::new(
        &transcript_models_map,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        10,
        8,
        5,
        2,
        5,
        20,
        1_000,
        20
    );

    assert_eq!(filter.identify_template_switch(&variant_call), Some(TemplateSwitchReason::StrongHomology));
    assert_eq!(filter.passes(&variant_call), false);
}

/// The same junction called a fusion gene: flagged, but only the pooled breakend family is
/// suppressed, so it passes.
#[test]
fn template_switch_filter_keeps_a_fusion_gene_with_strong_homology() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let gtf_file: String = directory.path().join("annotation.gtf").to_str().unwrap().to_string();
    // One gene, on a contig away from the junction: no exon edge at either breakend.
    std::fs::write(&gtf_file, concat!(
        "chrZ\tTEST\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_type \"protein_coding\"; level 1;\n",
        "chrZ\tTEST\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; transcript_type \"protein_coding\"; level 1;\n",
        "chrZ\tTEST\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_id \"E1\"; exon_number 1; level 1;\n"
    )).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(&gtf_file, "test", "v0");
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);
    let transcript_models_map: HashMap<ReadID, Arc<TranscriptModel>> = HashMap::new();

    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for read_id in [1, 2] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, 30, Strand::Forward, GraphOperationType::Downstream,
                1, 51, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::FusionGene
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);
    let filter: TemplateSwitchFilter<Gencode> = TemplateSwitchFilter::new(
        &transcript_models_map,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        10,
        8,
        5,
        2,
        5,
        20,
        1_000,
        20
    );

    assert_eq!(filter.identify_template_switch(&variant_call), None);
    assert_eq!(filter.passes(&variant_call), true);
}
