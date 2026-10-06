use bimap::BiMap;
use exacto_core::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

use super::*;


#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_tandem_duplication_1() {
    let chromosome_id: u16 = 0;
    let duplication_start: u32 = 7_673_001;
    let duplication_end: u32 = 7_673_500;
    let insertion_left_anchor: u32 = duplication_end;
    let insertion_right_anchor: u32 = insertion_left_anchor + 1;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let duplicated_sequence = fasta_map
        .get_sequence(
            "chr17",
            duplication_start as usize,
            duplication_end as usize,
        )
        .to_ascii_uppercase();

    assert_eq!(duplicated_sequence.len(), 500);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            duplicated_sequence.into_boxed_str(),
            VariantType::Insertion,
        ),
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![duplication_start, duplication_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records: Vec<VariantRecord> = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("Both ends of the duplicated insertion should align.");

    assert_eq!(retyped_records.len(), 1);

    let record = &retyped_records[0];

    assert!(record.get_position_1().abs_diff(duplication_start) < 10);
    assert_eq!(record.get_operation_1(), &GraphOperationType::Upstream);

    assert!(record.get_position_2().abs_diff(duplication_end) < 10);
    assert_eq!(record.get_operation_2(), &GraphOperationType::Downstream);

    assert_eq!(record.get_variant_type(), &VariantType::Breakpoint);
    assert!(record.get_sequence().is_empty());
}

#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_tandem_duplication_2() {
    let chromosome_id: u16 = 0;
    let duplication_start: u32 = 7_673_001;
    let duplication_end: u32 = 7_673_550;
    let insertion_left_anchor: u32 = duplication_end - 50; // insertion is deliberately misplaced
    let insertion_right_anchor: u32 = insertion_left_anchor + 1;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let duplicated_sequence = fasta_map
        .get_sequence(
            "chr17",
            duplication_start as usize,
            duplication_end as usize,
        )
        .to_ascii_uppercase();

    assert_eq!(duplicated_sequence.len(), 550);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            duplicated_sequence.into_boxed_str(),
            VariantType::Insertion,
        ),
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![duplication_start, duplication_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records: Vec<VariantRecord> = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("Both ends of the duplicated insertion should align.");

    assert_eq!(retyped_records.len(), 1);

    let record = &retyped_records[0];

    assert!(record.get_position_1().abs_diff(duplication_start) < 10);
    assert_eq!(record.get_operation_1(), &GraphOperationType::Upstream);

    assert!(record.get_position_2().abs_diff(duplication_end) < 10);
    assert_eq!(record.get_operation_2(), &GraphOperationType::Downstream);

    assert_eq!(record.get_variant_type(), &VariantType::Breakpoint);
    assert!(record.get_sequence().is_empty());
}

#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_inverted_duplication_1() {
    let chromosome_id: u16 = 0;
    let duplication_start: u32 = 7_673_001;
    let duplication_end: u32 = 7_673_500;
    let insertion_left_anchor: u32 = duplication_end;
    let insertion_right_anchor: u32 = insertion_left_anchor + 1;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let mut duplicated_sequence: String = fasta_map
        .get_sequence(
            "chr17",
            duplication_start as usize,
            duplication_end as usize,
        )
        .to_ascii_uppercase();
    duplicated_sequence = reverse_complement(&duplicated_sequence).to_string();

    assert_eq!(duplicated_sequence.len(), 500);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            duplicated_sequence.into_boxed_str(),
            VariantType::Insertion,
        ),
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![duplication_start, duplication_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records: Vec<VariantRecord> = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("Both ends of the duplicated insertion should align.");

    // We expect to get 2 breakpoints and one of them should be a fold-back junction.
    assert_eq!(retyped_records.len(), 2);

    // Junction 1: the fold-back junction. Base `duplication_end` is read twice, once per strand:
    // the forward pass leaves its downstream side into the turn, the reverse pass
    // enters it from the same side. Same position, both Downstream; the opposing
    // strands are what encode the inversion.
    let fold: &VariantRecord = &retyped_records[0];
    assert_eq!(*fold.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(fold.get_chromosome_1(), chromosome_id);
    assert_eq!(fold.get_position_1(), duplication_end);
    assert_eq!(*fold.get_strand_1(), Strand::Forward);
    assert_eq!(*fold.get_operation_1(), GraphOperationType::Downstream);
    assert_eq!(fold.get_chromosome_2(), chromosome_id);
    assert_eq!(fold.get_position_2(), duplication_end);
    assert_eq!(*fold.get_strand_2(), Strand::Reverse);
    assert_eq!(*fold.get_operation_2(), GraphOperationType::Downstream);
    assert_eq!(fold.get_sequence(), "");
    assert_eq!(fold.get_read_id(), 1);
    assert_eq!(fold.get_read_position_1(), 1_000);
    assert_eq!(fold.get_read_position_2(), 1_001);

    // Junction 2: the far junction. The inverted copy's other end (the segment's
    // LOW edge, traversed on the reverse strand, junction opening upstream) rejoins
    // the reference at the right insertion anchor.
    let far: &VariantRecord = &retyped_records[1];
    assert_eq!(*far.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(far.get_chromosome_1(), chromosome_id);
    assert_eq!(far.get_position_1(), duplication_start);
    assert_eq!(*far.get_strand_1(), Strand::Reverse);
    assert_eq!(*far.get_operation_1(), GraphOperationType::Upstream);
    assert_eq!(far.get_chromosome_2(), chromosome_id);
    assert_eq!(far.get_position_2(), insertion_right_anchor);
    assert_eq!(*far.get_strand_2(), Strand::Forward);
    assert_eq!(*far.get_operation_2(), GraphOperationType::Upstream);
    assert_eq!(far.get_sequence(), "");
    assert_eq!(far.get_read_id(), 1);
    assert_eq!(far.get_read_position_1(), 1_499);
    assert_eq!(far.get_read_position_2(), 1_500);
}

#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_inverted_duplication_2() {
    let chromosome_id: u16 = 0;
    let duplication_start: u32 = 7_673_001;
    let duplication_end: u32 = 7_673_550;
    let insertion_left_anchor: u32 = duplication_end - 50; // insertion is deliberately misplaced
    let insertion_right_anchor: u32 = insertion_left_anchor + 1;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let mut duplicated_sequence: String = fasta_map
        .get_sequence(
            "chr17",
            duplication_start as usize,
            duplication_end as usize,
        )
        .to_ascii_uppercase();
    duplicated_sequence = reverse_complement(&duplicated_sequence).to_string();

    assert_eq!(duplicated_sequence.len(), 550);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            duplicated_sequence.into_boxed_str(),
            VariantType::Insertion,
        ),
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![duplication_start, duplication_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records: Vec<VariantRecord> = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("Both ends of the duplicated insertion should align.");

    // We expect to get 2 breakpoints and one of them should be a fold-back junction.
    assert_eq!(retyped_records.len(), 2);

    // Junction 1: the fold-back junction. Base `duplication_end` is read twice, once per strand:
    // the forward pass leaves its downstream side into the turn, the reverse pass
    // enters it from the same side. Same position, both Downstream; the opposing
    // strands are what encode the inversion.
    let fold: &VariantRecord = &retyped_records[0];
    assert_eq!(*fold.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(fold.get_chromosome_1(), chromosome_id);
    assert!(fold.get_position_1().abs_diff(duplication_end) < 100);
    assert_eq!(*fold.get_strand_1(), Strand::Forward);
    assert_eq!(*fold.get_operation_1(), GraphOperationType::Downstream);
    assert_eq!(fold.get_chromosome_2(), chromosome_id);
    assert!(fold.get_position_2().abs_diff(duplication_end) < 100);
    assert_eq!(*fold.get_strand_2(), Strand::Reverse);
    assert_eq!(*fold.get_operation_2(), GraphOperationType::Downstream);
    assert_eq!(fold.get_sequence(), "");
    assert_eq!(fold.get_read_id(), 1);
    assert_eq!(fold.get_read_position_1(), 1_000);
    assert_eq!(fold.get_read_position_2(), 1_001);

    // Junction 2: the far junction. The inverted copy's other end (the segment's
    // LOW edge, traversed on the reverse strand, junction opening upstream) rejoins
    // the reference at the right insertion anchor.
    let far: &VariantRecord = &retyped_records[1];
    assert_eq!(*far.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(far.get_chromosome_1(), chromosome_id);
    assert!(far.get_position_1().abs_diff(duplication_start) < 100);
    assert_eq!(*far.get_strand_1(), Strand::Reverse);
    assert_eq!(*far.get_operation_1(), GraphOperationType::Upstream);
    assert_eq!(far.get_chromosome_2(), chromosome_id);
    assert!(far.get_position_2().abs_diff(insertion_right_anchor) < 100);
    assert_eq!(*far.get_strand_2(), Strand::Forward);
    assert_eq!(*far.get_operation_2(), GraphOperationType::Upstream);
    assert_eq!(far.get_sequence(), "");
    assert_eq!(far.get_read_id(), 1);
    assert_eq!(far.get_read_position_1(), 1_499);
    assert_eq!(far.get_read_position_2(), 1_500);
}

#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_inversion_1() {
    let chromosome_id: u16 = 0;
    let inversion_start: u32 = 7_673_001;
    let inversion_end: u32 = 7_673_500;

    let insertion_left_anchor: u32 = inversion_start - 1;   // 7,673,000
    let insertion_right_anchor: u32 = inversion_end + 1;    // 7,673,501

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let inverted_sequence: String = reverse_complement(
        fasta_map
            .get_sequence(
                "chr17",
                inversion_start as usize,
                inversion_end as usize,
            )
            .to_ascii_uppercase()
            .as_str()
    ).to_string();

    assert_eq!(inverted_sequence.len(), 500);

    let insertion_record: VariantRecord = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            inverted_sequence.into_boxed_str(),
            VariantType::Insertion
        )
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![inversion_start, inversion_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("Both ends of the inverted segment should align.");

    assert_eq!(retyped_records.len(), 2);

    // Junction 1: the low-side inversion breakend: the reference leaves the last
    // retained base heading downstream and enters the inverted segment at its HIGH
    // edge, also downstream. (L-1 D) join (H D).
    let low_side: &VariantRecord = &retyped_records[0];
    assert_eq!(*low_side.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(low_side.get_chromosome_1(), chromosome_id);
    assert_eq!(low_side.get_position_1(), insertion_left_anchor);
    assert_eq!(*low_side.get_strand_1(), Strand::Forward);
    assert_eq!(*low_side.get_operation_1(), GraphOperationType::Downstream);
    assert_eq!(low_side.get_chromosome_2(), chromosome_id);
    assert_eq!(low_side.get_position_2(), inversion_end);
    assert_eq!(*low_side.get_strand_2(), Strand::Reverse);
    assert_eq!(*low_side.get_operation_2(), GraphOperationType::Downstream);
    assert_eq!(low_side.get_sequence(), "");
    assert_eq!(low_side.get_read_position_1(), 1_000);
    assert_eq!(low_side.get_read_position_2(), 1_001);

    // Junction 2: the high-side inversion breakend. The read exits the inverted
    // segment at its LOW edge heading upstream and rejoins the reference at the first
    // retained base past the span. (L U) join (H+1 U).
    let high_side: &VariantRecord = &retyped_records[1];
    assert_eq!(*high_side.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(high_side.get_chromosome_1(), chromosome_id);
    assert_eq!(high_side.get_position_1(), inversion_start);
    assert_eq!(*high_side.get_strand_1(), Strand::Reverse);
    assert_eq!(*high_side.get_operation_1(), GraphOperationType::Upstream);
    assert_eq!(high_side.get_chromosome_2(), chromosome_id);
    assert_eq!(high_side.get_position_2(), insertion_right_anchor);
    assert_eq!(*high_side.get_strand_2(), Strand::Forward);
    assert_eq!(*high_side.get_operation_2(), GraphOperationType::Upstream);
    assert_eq!(high_side.get_sequence(), "");
    assert_eq!(high_side.get_read_position_1(), 1_499);
    assert_eq!(high_side.get_read_position_2(), 1_500);
}

#[test]
fn retype_insertion_as_breakpoints_returns_two_records_for_inversion_2() {
    let chromosome_id: u16 = 0;
    let inversion_start: u32 = 7_673_001;
    let inversion_end: u32 = 7_673_500;

    let insertion_left_anchor: u32 = inversion_start - 1 - 50;   // 7,673,000 - 50 (deliberately misplaced)
    let insertion_right_anchor: u32 = inversion_end + 1;         // 7,673,501

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let inverted_sequence: String = reverse_complement(
        fasta_map
            .get_sequence(
                "chr17",
                inversion_start as usize,
                inversion_end as usize,
            )
            .to_ascii_uppercase()
            .as_str()
    ).to_string();

    assert_eq!(inverted_sequence.len(), 500);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            inverted_sequence.into_boxed_str(),
            VariantType::Insertion
        )
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![inversion_start, inversion_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    )
        .expect("Both ends of the inverted segment should align.");

    assert_eq!(retyped_records.len(), 2);

    // Junction 1: the low-side inversion breakend: the reference leaves the last
    // retained base heading downstream and enters the inverted segment at its HIGH
    // edge, also downstream. (L-1 D) join (H D).
    let low_side: &VariantRecord = &retyped_records[0];
    assert_eq!(*low_side.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(low_side.get_chromosome_1(), chromosome_id);
    assert!(low_side.get_position_1().abs_diff(insertion_left_anchor) < 100);
    assert_eq!(*low_side.get_strand_1(), Strand::Forward);
    assert_eq!(*low_side.get_operation_1(), GraphOperationType::Downstream);
    assert_eq!(low_side.get_chromosome_2(), chromosome_id);
    assert!(low_side.get_position_2().abs_diff(inversion_end) < 100);
    assert_eq!(*low_side.get_strand_2(), Strand::Reverse);
    assert_eq!(*low_side.get_operation_2(), GraphOperationType::Downstream);
    assert_eq!(low_side.get_sequence(), "");
    assert_eq!(low_side.get_read_position_1(), 1_000);
    assert_eq!(low_side.get_read_position_2(), 1_001);

    // Junction 2: the high-side inversion breakend: the read exits the inverted
    // segment at its LOW edge heading upstream and rejoins the reference at the first
    // retained base past the span. (L U) join (H+1 U).
    let high_side: &VariantRecord = &retyped_records[1];
    assert_eq!(*high_side.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(high_side.get_chromosome_1(), chromosome_id);
    assert!(high_side.get_position_1().abs_diff(inversion_start) < 100);
    assert_eq!(*high_side.get_strand_1(), Strand::Reverse);
    assert_eq!(*high_side.get_operation_1(), GraphOperationType::Upstream);
    assert_eq!(high_side.get_chromosome_2(), chromosome_id);
    assert!(high_side.get_position_2().abs_diff(insertion_right_anchor) < 100);
    assert_eq!(*high_side.get_strand_2(), Strand::Forward);
    assert_eq!(*high_side.get_operation_2(), GraphOperationType::Upstream);
    assert_eq!(high_side.get_sequence(), "");
    assert_eq!(high_side.get_read_position_1(), 1_499);
    assert_eq!(high_side.get_read_position_2(), 1_500);
}

#[test]
fn retype_insertion_as_breakpoints_returns_one_record_for_deletion() {
    let chromosome_id: u16 = 0;
    let insertion_left_anchor: u32 = 2_000_000;
    let insertion_right_anchor: u32 = 2_000_001;
    let donor_start: u32 = 3_500_001;
    let donor_end: u32 = 3_500_500;

    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map = FastaMap::new(fasta_full_path.to_str().unwrap());

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), chromosome_id);

    let donor_sequence: String = fasta_map
        .get_sequence("chr17", donor_start as usize, donor_end as usize)
        .to_ascii_uppercase();

    assert_eq!(donor_sequence.len(), 500);

    let insertion_record = VariantRecord::new(
        1,
        1_000,
        1_500,
        GraphOperation::new(
            chromosome_id,
            insertion_left_anchor,
            Strand::Forward,
            GraphOperationType::Downstream,
            chromosome_id,
            insertion_right_anchor,
            Strand::Forward,
            GraphOperationType::Upstream,
            donor_sequence.into_boxed_str(),
            VariantType::Insertion
        )
    );

    let query_positions: HashMap<u16, Vec<u32>> = HashMap::from([(
        chromosome_id,
        vec![donor_start, donor_end]
    )]);

    let options: IdentifyGermlineDNAVariantsOptions = IdentifyGermlineDNAVariantsOptions::default();

    let retyped_records = retype_insertion_as_breakpoints(
        &insertion_record,
        &query_positions,
        &fasta_map,
        &chromosome_names_map,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_max_ins_len,
        options.calling.bkpt_rescue_search_distance,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_span_proportion
    ).expect("The far junction should survive.");

    assert_eq!(retyped_records.len(), 1);

    let kept: &VariantRecord = &retyped_records[0];
    assert_eq!(*kept.get_variant_type(), VariantType::Breakpoint);
    assert_eq!(kept.get_chromosome_1(), chromosome_id);
    assert_eq!(kept.get_position_1(), insertion_right_anchor);
    assert_eq!(*kept.get_strand_1(), Strand::Forward);
    assert_eq!(*kept.get_operation_1(), GraphOperationType::Upstream);
    assert_eq!(kept.get_chromosome_2(), chromosome_id);
    assert_eq!(kept.get_position_2(), donor_end);
    assert_eq!(*kept.get_strand_2(), Strand::Forward);
    assert_eq!(*kept.get_operation_2(), GraphOperationType::Downstream);
    assert_eq!(kept.get_sequence(), "");
    assert_eq!(kept.get_read_position_1(), 1_499);
    assert_eq!(kept.get_read_position_2(), 1_500);

    // Explicit negative: nothing may spell the phantom deletion's breakends.
    assert!(retyped_records.iter().all(|record| {
        record.get_position_1() != insertion_left_anchor
            && record.get_position_2() != donor_start
    }));
}
#[test]
fn place_query_returns_the_best_placement_whatever_the_order_of_the_windows() {
    // One contig holding a 200-base sequence at 500-699 and a copy of it with 17 substitutions
    // at 2000-2199. The query is the sequence.
    let mut state: u32 = 2_463_534_242;
    let mut contig: Vec<u8> = (0..3_000)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state % 4) as usize]
        })
        .collect();
    let query: Vec<u8> = contig[499..699].to_vec();
    for i in 0..200 {
        contig[1_999 + i] = if i % 12 == 5 {
            match query[i] { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' }
        } else {
            query[i]
        };
    }
    let temp_dir = tempfile::tempdir().unwrap();
    let fasta_file = temp_dir.path().join("two_copies.fa");
    let lines: Vec<String> = contig.chunks(60).map(|line| String::from_utf8(line.to_vec()).unwrap()).collect();
    fs::write(&fasta_file, format!(">syn\n{}\n", lines.join("\n"))).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("syn".into(), 0);

    let options: DNAVariantCallingOptions = DNAVariantCallingOptions::DEFAULT;
    let query: String = String::from_utf8(query).unwrap();
    let place = |windows: &[(u16, u32, u32)]| -> Option<(u32, u32, Strand, i32)> {
        place_query(
            &query,
            windows,
            &fasta_map,
            &chromosome_names_map,
            options.bkpt_rescue_realignment_gap_open_score,
            options.bkpt_rescue_realignment_gap_extend_score,
            options.bkpt_rescue_realignment_k,
            options.bkpt_rescue_realignment_band_width,
            options.bkpt_rescue_realignment_min_score_fraction,
            options.bkpt_rescue_realignment_min_query_coverage
        ).map(|placement| (placement.low, placement.high, placement.strand, placement.score))
    };
    let exact: (u16, u32, u32) = (0, 400, 800);
    let copy: (u16, u32, u32) = (0, 1_900, 2_300);

    assert_eq!(place(&[copy]), Some((2_000, 2_199, Strand::Forward, 166)));
    assert_eq!(place(&[exact, copy]), Some((500, 699, Strand::Forward, 200)));
    assert_eq!(place(&[copy, exact]), Some((500, 699, Strand::Forward, 200)));
}

#[test]
fn place_query_returns_none_for_two_placements_of_equal_score() {
    // One contig holding the same 200-base sequence at 500-699 and at 2000-2199.
    let mut state: u32 = 2_463_534_242;
    let mut contig: Vec<u8> = (0..3_000)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state % 4) as usize]
        })
        .collect();
    let query: Vec<u8> = contig[499..699].to_vec();
    contig[1_999..2_199].copy_from_slice(&query);
    let temp_dir = tempfile::tempdir().unwrap();
    let fasta_file = temp_dir.path().join("two_copies.fa");
    let lines: Vec<String> = contig.chunks(60).map(|line| String::from_utf8(line.to_vec()).unwrap()).collect();
    fs::write(&fasta_file, format!(">syn\n{}\n", lines.join("\n"))).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("syn".into(), 0);

    let options: DNAVariantCallingOptions = DNAVariantCallingOptions::DEFAULT;
    let query: String = String::from_utf8(query).unwrap();
    let place = |windows: &[(u16, u32, u32)]| -> Option<(u32, u32)> {
        place_query(
            &query,
            windows,
            &fasta_map,
            &chromosome_names_map,
            options.bkpt_rescue_realignment_gap_open_score,
            options.bkpt_rescue_realignment_gap_extend_score,
            options.bkpt_rescue_realignment_k,
            options.bkpt_rescue_realignment_band_width,
            options.bkpt_rescue_realignment_min_score_fraction,
            options.bkpt_rescue_realignment_min_query_coverage
        ).map(|placement| (placement.low, placement.high))
    };
    let first: (u16, u32, u32) = (0, 400, 800);
    let second: (u16, u32, u32) = (0, 1_900, 2_300);

    assert_eq!(place(&[first]), Some((500, 699)));
    assert_eq!(place(&[first, second]), None);
    assert_eq!(place(&[second, first]), None);
}
