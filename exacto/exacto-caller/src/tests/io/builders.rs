use bimap::BiMap;
use exacto_core::prelude::Strand;
use std::collections::HashSet;

use super::*;


#[test]
fn build_dna_variant_records_leaves_the_size_of_a_translocation_empty() {
    let translocation: GraphOperation = GraphOperation::new(
        0,
        7_676_155,
        Strand::Forward,
        GraphOperationType::Downstream,
        1,
        5_170_101,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Translocation
    );
    let snv: GraphOperation = GraphOperation::new(
        0,
        7_674_224,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        7_674_226,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::SingleNucleotideVariant
    );
    let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new(DNAVariantOrigin::Somatic);
    for (id, operation) in [(1, translocation), (2, snv)] {
        let variant_records: HashSet<VariantRecord> = HashSet::from([VariantRecord::new(1, 99, 100, operation)]);
        let mut variant_call: VariantCall = VariantCall::from_variant_records(id, variant_records, 1, -1, -1, -1);
        variant_call.set_id(id);
        variant_call_set.add_variant_call(variant_call);
    }
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0);
    chromosome_names_map.insert("chr18".into(), 1);
    variant_call_set.load_chromosome_names(chromosome_names_map);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    read_names_map.insert("read-1".into(), 1);
    variant_call_set.load_read_names(read_names_map);

    let sizes: Vec<(Box<str>, Option<i32>)> = build_dna_variant_records(&variant_call_set)
        .map(|record| (record.variant_type, record.variant_size))
        .collect();

    assert_eq!(sizes, vec![("TRA".into(), None), ("SNV".into(), Some(1))]);
}


/// The scga-mini DNA and RNA samples carry the same variants: rna-001's SNV is dna-001's, and
/// the germline (tumor alone) and somatic tables of dna-001 both hold it. rna-005's non-canonical
/// splice is dna-005's deletion at the same breakpoints, but the types differ: no match.
#[test]
fn annotate_dna_variant_origins_matches_the_scga_mini_dna_calls() {
    let data_dir = std::path::Path::new(env!("EXACTO_TEST_DATA"));
    for (sample, expected) in [("001", vec![("SNV", "germline;somatic")]), ("005", vec![("NCS", "")])] {
        let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
        for origin in ["germline", "somatic"] {
            dna_variant_records.extend(load_dna_variant_records(
                data_dir.join(format!("variant_calling/dna/scga-mini-dna-{}-tumor_exacto_{}_dna_variants.tsv", sample, origin)).to_str().unwrap()
            ));
        }
        let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = load_assembled_transcript_variant_records(
            data_dir.join(format!("variant_calling/rna/scga-mini-rna-{}-tumor_exacto_assembled_transcript_variants.tsv", sample)).to_str().unwrap()
        );
        let origins: Vec<(Box<str>, Box<str>)> = annotate_dna_variant_origins(rna_variant_records, &dna_variant_records, 1)
            .map(|record| (record.variant_type, record.origin))
            .collect();
        assert_eq!(
            origins,
            expected.into_iter().map(|(variant_type, origin)| (variant_type.into(), origin.into())).collect::<Vec<(Box<str>, Box<str>)>>(),
            "scga-mini-rna-{}", sample
        );
    }
}

/// A match is the same chromosomes, operations, type and sequence (any case), each position
/// within the buffer; the strand is not compared.
#[test]
fn annotate_dna_variant_origins_compares_everything_but_the_strand() {
    let dna_variant_records: Vec<DNAVariantRecord> = vec![
        DNAVariantRecord {
            variant_id: 1,
            origin: "somatic".into(),
            chromosome_1: "chr17".into(),
            position_1: 7_674_224,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7_674_226,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            ..Default::default()
        }
    ];
    // (position_1, position_2, variant_type, sequence, buffer, expected origin)
    let cases: Vec<(u32, u32, &str, &str, u32, &str)> = vec![
        (7_674_224, 7_674_226, "SNV", "a", 0, "somatic"),
        (7_674_225, 7_674_227, "SNV", "A", 0, ""),
        (7_674_225, 7_674_227, "SNV", "A", 1, "somatic"),
        (7_674_223, 7_674_225, "SNV", "A", 1, "somatic"),
        (7_674_226, 7_674_228, "SNV", "A", 1, ""),
        (7_674_224, 7_674_226, "SNV", "C", 1, ""),
        (7_674_224, 7_674_226, "MNV", "A", 1, "")
    ];
    for (position_1, position_2, variant_type, sequence, buffer, expected) in cases {
        let rna_variant_record: AssembledTranscriptVariantRecord = AssembledTranscriptVariantRecord {
            variant_id: 1,
            assembled_transcript_name: "transcript".into(),
            chromosome_1: "chr17".into(),
            position_1: position_1,
            strand_1: "-".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: position_2,
            strand_2: "-".into(),
            operation_2: "U".into(),
            variant_size: Some(1),
            variant_type: variant_type.into(),
            sequence: sequence.into(),
            ..Default::default()
        };
        let origins: Vec<Box<str>> = annotate_dna_variant_origins(vec![rna_variant_record], &dna_variant_records, buffer)
            .map(|record| record.origin)
            .collect();
        assert_eq!(origins, vec![Box::<str>::from(expected)], "{} {} {} {} buffer {}", position_1, position_2, variant_type, sequence, buffer);
    }

    // No DNA variants: every origin is empty.
    let rna_variant_record: AssembledTranscriptVariantRecord = AssembledTranscriptVariantRecord {
        chromosome_1: "chr17".into(),
        position_1: 7_674_224,
        chromosome_2: "chr17".into(),
        position_2: 7_674_226,
        variant_type: "SNV".into(),
        sequence: "A".into(),
        ..Default::default()
    };
    let origins: Vec<Box<str>> = annotate_dna_variant_origins(vec![rna_variant_record], &[], 1)
        .map(|record| record.origin)
        .collect();
    assert_eq!(origins, vec![Box::<str>::from("")]);
}
