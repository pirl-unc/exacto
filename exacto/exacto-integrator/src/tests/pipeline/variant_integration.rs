use std::collections::HashSet;
use exacto_caller::prelude::*;
use std::fs;
use std::path::Path;
use exacto_core::prelude::Gencode;

use super::*;


#[test]
fn test_variant_integration_1() {
    // scga-mini case 001: the somatic TP53 SNV at chr17:7,674,225, called in DNA and RNA.
    // The expected rows are the integrate-vars output of scripts/data/06_integration.
    let tsv_path_1 = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-001-tumor_exacto_somatic_dna_variants.tsv");
    let tsv_full_path_1 = fs::canonicalize(tsv_path_1).unwrap();
    let tsv_file_1: &str = tsv_full_path_1.to_str().unwrap();
    let tsv_path_2 = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/rna/scga-mini-rna-001-tumor_exacto_assembled_transcript_variants.tsv");
    let tsv_full_path_2 = fs::canonicalize(tsv_path_2).unwrap();
    let tsv_file_2: &str = tsv_full_path_2.to_str().unwrap();
    let tsv_path_3 = Path::new(env!("EXACTO_TEST_DATA")).join("integration/scga-mini-rna-001-tumor_exacto_integrated_variants.tsv");
    let tsv_full_path_3 = fs::canonicalize(tsv_path_3).unwrap();
    let tsv_file_3: &str = tsv_full_path_3.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();
    let dna_variant_records: Vec<DNAVariantRecord> = load_dna_variant_records(tsv_file_1);
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = load_assembled_transcript_variant_records(tsv_file_2);
    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );
    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        100_000,
        1
    ).unwrap();

    assert_eq!(integrated_variants.len(), 1);
    assert_eq!(integrated_variants[0].rna_variant_id, 1);
    assert_eq!(integrated_variants[0].dna_variant_ids.keys().copied().collect::<Vec<u32>>(), vec![1]);
    assert_eq!(integrated_variants[0].dna_variant_ids[&1].get_distance(), 0);

    let integrated_variant_records: Vec<IntegratedVariantRecord> = build_integrated_variant_records(&integrated_variants).collect();
    assert_eq!(integrated_variant_records, load_integrated_variant_records(tsv_file_3));
}

#[test]
fn test_variant_integration_2() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    let dna_variant_record = DNAVariantRecord {
        origin: "somatic".into(),
        variant_id: 1,
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
        consensus_read_names: "read1".into(),
        num_consensus_read_names: 1,
        read_names: "read1".into(),
        num_read_names: 1
    };
    dna_variant_records.push(dna_variant_record);

    let mut rna_variant_records: Vec<AssembledTranscriptVariantRecord> = Vec::new();
    let rna_variant_record = AssembledTranscriptVariantRecord {
        variant_id: 1,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "TP53".into(),
        reference_transcript_id: "ENST00000269305.9".into(),
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
        read_start: 100,
        read_end: 100,
        origin: "".into()
    };
    rna_variant_records.push(rna_variant_record);

    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        1000,
        100_000,
        1
    ).unwrap();

    assert_eq!(integrated_variants.len(), 1);
}


#[test]
fn minus_strand_rna_variant_past_3prime_end_is_linked_to_last_exon() {
    // TP53 is on the minus strand: exon 11 is at transcript.start (7,668,421), exon 1 at
    // transcript.end (7,687,490). The RNA SNV lies 2 kb below transcript.start.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    for (variant_id, position) in [(1, 7_668_900), (2, 7_687_400)] {
        dna_variant_records.push(DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: variant_id,
            chromosome_1: "chr17".into(),
            position_1: position - 1,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: position + 1,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read1".into(),
            num_consensus_read_names: 1,
            read_names: "read1".into(),
            num_read_names: 1
        });
    }
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "TP53".into(),
        reference_transcript_id: "ENST00000269305.9".into(),
        chromosome_1: "chr17".into(),
        position_1: 7_666_420,
        strand_1: "-".into(),
        operation_1: "D".into(),
        chromosome_2: "chr17".into(),
        position_2: 7_666_422,
        strand_2: "-".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        read_start: 100,
        read_end: 100,
        origin: "".into()
    }];

    // A 1 kb intergenic distance keeps the distance rule out of the way.
    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        1_000,
        1
    ).unwrap();

    assert_eq!(integrated_variants.len(), 1);
    let dna_variant_ids: Vec<u32> = integrated_variants[0].dna_variant_ids.keys().copied().collect();
    assert_eq!(dna_variant_ids, vec![1]);
    assert_eq!(integrated_variants[0].dna_variant_ids[&1].get_distance(), 2_477);
}

#[test]
fn rna_variant_moving_away_from_transcript_never_gains_a_link() {
    // RNA SNVs above the 5' end of TP53 (exon 1 ends at 7,687,490). DNA SNV 1 is in exon 5,
    // DNA SNV 2 in exon 2, and DNA SNV 3 is intergenic, 30 kb above the 5' end.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    for (variant_id, position) in [(1, 7_675_100), (2, 7_676_550), (3, 7_717_490)] {
        dna_variant_records.push(DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: variant_id,
            chromosome_1: "chr17".into(),
            position_1: position - 1,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: position + 1,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read1".into(),
            num_consensus_read_names: 1,
            read_names: "read1".into(),
            num_read_names: 1
        });
    }
    let offsets: [u32; 7] = [500, 5_000, 9_000, 11_000, 15_000, 40_000, 59_000];
    let mut rna_variant_records: Vec<AssembledTranscriptVariantRecord> = Vec::new();
    for (i, offset) in offsets.iter().enumerate() {
        rna_variant_records.push(AssembledTranscriptVariantRecord {
            variant_id: i as u32 + 1,
            assembled_transcript_name: "transcript_1".into(),
            reference_gene_name: "TP53".into(),
            reference_transcript_id: "ENST00000269305.9".into(),
            chromosome_1: "chr17".into(),
            position_1: 7_687_490 + offset - 1,
            strand_1: "-".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7_687_490 + offset + 1,
            strand_2: "-".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            read_start: 100,
            read_end: 100,
            origin: "".into()
        });
    }

    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        100_000,
        1
    ).unwrap();

    // Exon 5 is four exons from exon 1: never linked. Exon 2 is linked while the RNA SNV is
    // within 10 kb of the 5' end. The intergenic SNV is within 100 kb of every RNA SNV.
    assert_eq!(integrated_variants.len(), offsets.len());
    for (integrated_variant, offset) in integrated_variants.iter().zip(offsets.iter()) {
        let dna_variant_ids: Vec<u32> = integrated_variant.dna_variant_ids.keys().copied().collect();
        let expected: Vec<u32> = if *offset < 10_000 { vec![2, 3] } else { vec![3] };
        assert_eq!(dna_variant_ids, expected, "RNA SNV {} bp above TP53", offset);
    }
}

#[test]
fn rna_variant_on_another_chromosome_is_not_placed_on_transcript() {
    // The RNA SNV is on chr18 at the number of a base in TP53 exon 7, and lists a chr18
    // transcript (PTPRM) and TP53. The DNA SNV is at that base of chr17.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        origin: "somatic".into(),
        variant_id: 1,
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
        consensus_read_names: "read1".into(),
        num_consensus_read_names: 1,
        read_names: "read1".into(),
        num_read_names: 1
    }];
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "PTPRM;TP53".into(),
        reference_transcript_id: "ENST00000580170.6;ENST00000269305.9".into(),
        chromosome_1: "chr18".into(),
        position_1: 7_674_224,
        strand_1: "+".into(),
        operation_1: "D".into(),
        chromosome_2: "chr18".into(),
        position_2: 7_674_226,
        strand_2: "+".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        read_start: 100,
        read_end: 100,
        origin: "".into()
    }];

    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        100_000,
        1
    ).unwrap();

    assert_eq!(integrated_variants.len(), 0);
}

#[test]
fn fusion_is_not_linked_through_its_partner_chromosome() {
    // Side 1 of the fusion is on chr17, 180 kb above TP53. Side 2 is on chr18 at the number of
    // a base in TP53 exon 7. The DNA SNV is at that base of chr17.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        origin: "somatic".into(),
        variant_id: 1,
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
        consensus_read_names: "read1".into(),
        num_consensus_read_names: 1,
        read_names: "read1".into(),
        num_read_names: 1
    }];
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "TP53;PTPRM".into(),
        reference_transcript_id: "ENST00000269305.9;ENST00000580170.6".into(),
        chromosome_1: "chr17".into(),
        position_1: 7_867_490,
        strand_1: "-".into(),
        operation_1: "D".into(),
        chromosome_2: "chr18".into(),
        position_2: 7_674_225,
        strand_2: "+".into(),
        operation_2: "U".into(),
        sequence: "".into(),
        variant_size: None,
        variant_type: "FUS".into(),
        read_start: 100,
        read_end: 100,
        origin: "".into()
    }];

    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        100_000,
        1
    ).unwrap();

    assert_eq!(integrated_variants.len(), 0);
}

#[test]
fn terminal_exon_window_is_the_same_at_both_ends() {
    // One RNA SNV in each of the 11 exons of TP53. DNA SNV 1 is 1 kb above the 5' end
    // (exon 1), DNA SNV 2 is 1 kb below the 3' end (exon 11).
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    for (variant_id, position) in [(1, 7_688_490), (2, 7_667_421)] {
        dna_variant_records.push(DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: variant_id,
            chromosome_1: "chr17".into(),
            position_1: position - 1,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: position + 1,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read1".into(),
            num_consensus_read_names: 1,
            read_names: "read1".into(),
            num_read_names: 1
        });
    }
    let exon_middles: [u32; 11] = [
        7_687_433, 7_676_571, 7_676_392, 7_676_133, 7_675_144, 7_674_915,
        7_674_235, 7_673_769, 7_673_571, 7_670_662, 7_669_055
    ];
    let mut rna_variant_records: Vec<AssembledTranscriptVariantRecord> = Vec::new();
    for (i, position) in exon_middles.iter().enumerate() {
        rna_variant_records.push(AssembledTranscriptVariantRecord {
            variant_id: i as u32 + 1,
            assembled_transcript_name: "transcript_1".into(),
            reference_gene_name: "TP53".into(),
            reference_transcript_id: "ENST00000269305.9".into(),
            chromosome_1: "chr17".into(),
            position_1: position - 1,
            strand_1: "-".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: position + 1,
            strand_2: "-".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            read_start: 100,
            read_end: 100,
            origin: "".into()
        });
    }

    let integrated_variants: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records,
        &rna_variant_records,
        &gene_annotator,
        2,
        10_000,
        100_000,
        1
    ).unwrap();

    let links: Vec<(u32, Vec<u32>)> = integrated_variants
        .iter()
        .map(|integrated_variant| (integrated_variant.rna_variant_id, integrated_variant.dna_variant_ids.keys().copied().collect()))
        .collect();
    assert_eq!(links, vec![
        (1, vec![1]), (2, vec![1]), (3, vec![1]),
        (9, vec![2]), (10, vec![2]), (11, vec![2])
    ]);
}

#[test]
fn integration_lists_dna_variants_in_ascending_id_order() {
    // One RNA SNV in TP53 exon 7 and 12 DNA SNVs in the same exon, given in shuffled ID order.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    for (i, variant_id) in [7, 3, 11, 1, 9, 5, 12, 2, 8, 10, 4, 6].into_iter().enumerate() {
        dna_variant_records.push(DNAVariantRecord {
            origin: "somatic".into(),
            variant_id: variant_id,
            chromosome_1: "chr17".into(),
            position_1: 7_674_190 + 5 * i as u32,
            strand_1: "+".into(),
            operation_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2: 7_674_192 + 5 * i as u32,
            strand_2: "+".into(),
            operation_2: "U".into(),
            sequence: "A".into(),
            variant_size: Some(1),
            variant_type: "SNV".into(),
            consensus_read_names: "read1".into(),
            num_consensus_read_names: 1,
            read_names: "read1".into(),
            num_read_names: 1
        });
    }
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "TP53".into(),
        reference_transcript_id: "ENST00000269305.9".into(),
        chromosome_1: "chr17".into(),
        position_1: 7_674_224,
        strand_1: "-".into(),
        operation_1: "D".into(),
        chromosome_2: "chr17".into(),
        position_2: 7_674_226,
        strand_2: "-".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        read_start: 100,
        read_end: 100,
        origin: "".into()
    }];

    let integrated_variants_1: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records, &rna_variant_records, &gene_annotator, 2, 10_000, 100_000, 2
    ).unwrap();
    let integrated_variants_2: Vec<RNAVariantIntegration> = integrate_dna_rna_variants(
        &dna_variant_records, &rna_variant_records, &gene_annotator, 2, 10_000, 100_000, 2
    ).unwrap();

    let dna_variant_ids: Vec<u32> = build_integrated_variant_records(&integrated_variants_1)
        .map(|record| record.dna_variant_id)
        .collect();
    assert_eq!(dna_variant_ids, (1..=12).collect::<Vec<u32>>());
    let mut hasher_1 = std::collections::hash_map::DefaultHasher::new();
    let mut hasher_2 = std::collections::hash_map::DefaultHasher::new();
    std::hash::Hash::hash(&integrated_variants_1, &mut hasher_1);
    std::hash::Hash::hash(&integrated_variants_2, &mut hasher_2);
    assert_eq!(integrated_variants_1, integrated_variants_2);
    assert_eq!(std::hash::Hasher::finish(&hasher_1), std::hash::Hasher::finish(&hasher_2));
}

#[test]
fn unknown_reference_transcript_is_an_error() {
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 7,
        assembled_transcript_name: "transcript_1".into(),
        reference_gene_name: "TP53".into(),
        reference_transcript_id: "ENST00000269305.9;ENST00000000000.1".into(),
        chromosome_1: "chr17".into(),
        position_1: 7_674_224,
        strand_1: "-".into(),
        operation_1: "D".into(),
        chromosome_2: "chr17".into(),
        position_2: 7_674_226,
        strand_2: "-".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        read_start: 100,
        read_end: 100,
        origin: "".into()
    }];

    let result = integrate_dna_rna_variants(&Vec::new(), &rna_variant_records, &gene_annotator, 2, 10_000, 100_000, 1);

    match result {
        Err(IntegratorError::UnknownTranscript { transcript_id, rna_variant_id, .. }) => {
            assert_eq!(&*transcript_id, "ENST00000000000.1");
            assert_eq!(rna_variant_id, 7);
        },
        other => panic!("expected UnknownTranscript, got {:?}", other)
    }
}

#[test]
fn repeated_dna_variant_id_is_an_error() {
    // Germline and somatic calls each number their IDs from 1.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut dna_variant_records: Vec<DNAVariantRecord> = Vec::new();
    for origin in ["germline", "somatic"] {
        dna_variant_records.push(DNAVariantRecord {
            origin: origin.into(),
            variant_id: 1,
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
            consensus_read_names: "read1".into(),
            num_consensus_read_names: 1,
            read_names: "read1".into(),
            num_read_names: 1
        });
    }

    let result = integrate_dna_rna_variants(&dna_variant_records, &Vec::new(), &gene_annotator, 2, 10_000, 100_000, 1);

    match result {
        Err(IntegratorError::DuplicateDNAVariantId { variant_id, origin_1, origin_2 }) => {
            assert_eq!(variant_id, 1);
            assert_eq!((&*origin_1, &*origin_2), ("germline", "somatic"));
        },
        other => panic!("expected DuplicateDNAVariantId, got {:?}", other)
    }
}

#[test]
fn transcript_exons_locate_as_transcript_does() {
    // Every exon edge, the base on each side of it, and the bases beyond each transcript end,
    // for every transcript of the fixture.
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
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

    let mut num_positions: usize = 0;
    for transcript in gene_annotator.get_transcripts() {
        let transcript_exons = TranscriptExons::new(transcript);
        let mut positions: Vec<u32> = vec![transcript.start - 1, transcript.end + 1];
        for exon in transcript.exons.values() {
            positions.extend([exon.start - 1, exon.start, exon.end, exon.end + 1]);
        }
        for position in positions {
            let expected = match transcript.locate_position(position) {
                (GenicRegion::Exonic, Some((exon, _))) => format!("exonic {}", exon.exon_number),
                (GenicRegion::Intronic, Some((exon_1, exon_2))) => format!("intronic {} {}", exon_1.exon_number, exon_2.exon_number),
                _ => "intergenic".to_string()
            };
            let actual = match transcript_exons.locate(position) {
                Location::Exonic(exon_number) => format!("exonic {}", exon_number),
                Location::Intronic(exon_number_1, exon_number_2) => {
                    let (a, b) = if transcript.strand == Strand::Forward { (exon_number_1, exon_number_2) } else { (exon_number_2, exon_number_1) };
                    format!("intronic {} {}", a, b)
                },
                Location::Intergenic { .. } => "intergenic".to_string()
            };
            assert_eq!(actual, expected, "{} at {}", transcript.transcript_id, position);
            num_positions += 1;
        }
    }
    assert!(num_positions > 100_000);
}
