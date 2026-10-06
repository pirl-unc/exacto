use exacto_core::prelude::Strand;
use std::sync::Arc;

use super::*;



#[test]
fn diff_dna_variant_records_returns_matches_1() {
    // INS:chr1:1000
    // DEL:chr1:1100-1200
    // SNV:chr1:1205
    // SNV:chr1:1400
    let mut a: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1001
    // SNV:chr1:1205
    let mut b: Vec<VariantRecord> = Vec::new();

    // INS:chr1:1000
    let go_a1: GraphOperation = GraphOperation::new(
        0,
        1000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "ACGATCGACT".into(),
        VariantType::Insertion
    );
    let vr_a1: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_a1
    );

    // DEL:chr1:1101-1200
    let go_a2: GraphOperation = GraphOperation::new(
        0,
        1100,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1201,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_a2: VariantRecord = VariantRecord::new(
        1,
        2,
        104,
        go_a2
    );

    // SNV:chr1:1205
    let go_a3: GraphOperation = GraphOperation::new(
        0,
        1204,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_a3: VariantRecord = VariantRecord::new(
        1,
        105,
        107,
        go_a3
    );

    // SNV:chr1:1400
    let go_a4: GraphOperation = GraphOperation::new(
        0,
        1399,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1401,
        Strand::Forward,
        GraphOperationType::Upstream,
        "T".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_a4: VariantRecord = VariantRecord::new(
        1,
        108,
        110,
        go_a4
    );

    a.push(vr_a1);
    a.push(vr_a2);
    a.push(vr_a3);
    a.push(vr_a4);

    // INS:chr1:1001
    let go_b1: GraphOperation = GraphOperation::new(
        0,
        1001,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1002,
        Strand::Forward,
        GraphOperationType::Upstream,
        "CGATCGACTA".into(),
        VariantType::Insertion
    );
    let vr_b1: VariantRecord = VariantRecord::new(
        2,
        111,
        112,
        go_b1
    );

    // SNV:chr1:1205
    let go_b2: GraphOperation = GraphOperation::new(
        0,
        1204,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1206,
        Strand::Forward,
        GraphOperationType::Upstream,
        "A".into(),
        VariantType::SingleNucleotideVariant
    );
    let vr_b2: VariantRecord = VariantRecord::new(
        2,
        113,
        115,
        go_b2
    );

    b.push(vr_b1);
    b.push(vr_b2);

    let a_ref: Vec<Arc<VariantRecord>> = a
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let b_ref: Vec<Arc<VariantRecord>> = b
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let a_only: Vec<Arc<VariantRecord>> = diff_dna_variant_records(
        a_ref,
        b_ref,
        100_000,
        1,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        0,
        true,
        false
    );

    assert!(a_only.len() == 2);

    let mut found_del_1: bool = false;
    let mut found_snv_2: bool = false;

    for variant_record in a_only.iter() {
        if variant_record.get_graph_operation_boxed_str() == "0:1100:+:D:0:1201:+:U::0:DEL".into() {
            found_del_1 = true;
        }
        if variant_record.get_graph_operation_boxed_str() == "0:1399:+:D:0:1401:+:U:T:1:SNV".into() {
            found_snv_2 = true;
        }
    }

    assert!(found_del_1);
    assert!(found_snv_2);
}

#[test]
fn diff_dna_variant_records_returns_matches_2() {
    let mut a: Vec<VariantRecord> = Vec::new();
    let mut b: Vec<VariantRecord> = Vec::new();

    // DEL:chr1:1010001-1020000
    let go_a1: GraphOperation = GraphOperation::new(
        0,
        1010000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1020001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_a1: VariantRecord = VariantRecord::new(
        1,
        0,
        1,
        go_a1
    );
    a.push(vr_a1);

    // DEL:chr1:1005001-1020000
    let go_b1: GraphOperation = GraphOperation::new(
        0,
        1005000,
        Strand::Forward,
        GraphOperationType::Downstream,
        0,
        1020001,
        Strand::Forward,
        GraphOperationType::Upstream,
        "".into(),
        VariantType::Deletion
    );
    let vr_b1: VariantRecord = VariantRecord::new(
        2,
        2,
        3,
        go_b1
    );
    b.push(vr_b1);

    let a_ref: Vec<Arc<VariantRecord>> = a
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let b_ref: Vec<Arc<VariantRecord>> = b
        .iter()
        .map(|record| Arc::new(record.clone()))
        .collect();
    let a_only: Vec<Arc<VariantRecord>> = diff_dna_variant_records(
        a_ref,
        b_ref,
        100_000,
        1,
        0.5f64,
        0.5f64,
        1000,
        0.01f64,
        1,
        true,
        false
    );

    assert_eq!(a_only.is_empty(), true);
}


/// A tumour deletion whose position 2 lies within the clustering distance (1,000) of the contig
/// start, against a control that holds another deletion: the search around position 2 starts
/// at base 0 rather than below it.
#[test]
fn diff_dna_variant_records_searches_a_deletion_near_the_contig_start() {
    let deletion = |read_id: usize, position_1: u32, position_2: u32| -> Arc<VariantRecord> {
        Arc::new(VariantRecord::new(
            read_id,
            0,
            1,
            GraphOperation::new(
                0,
                position_1,
                Strand::Forward,
                GraphOperationType::Downstream,
                0,
                position_2,
                Strand::Forward,
                GraphOperationType::Upstream,
                "".into(),
                VariantType::Deletion
            )
        ))
    };
    let a: Vec<Arc<VariantRecord>> = vec![deletion(1, 300, 321)];
    let b: Vec<Arc<VariantRecord>> = vec![deletion(2, 5_000, 5_201)];

    let a_only: Vec<Arc<VariantRecord>> = diff_dna_variant_records(a, b, 100_000, 1, 0.5, 0.5, 1_000, 0.01, 0, true, false);

    assert_eq!(a_only.len(), 1);
    assert_eq!(a_only[0].get_position_2(), 321);
}


/// Two control reads hold the tumour's SNV, one read from each strand. They are two reads of one
/// allele, so a control record needs two reads (`max_control_reads` 1) and this one has them.
/// A record holds its bases in the orientation of its strand 1, so the reverse read spells the
/// alternate base A as T.
#[test]
fn diff_dna_variant_records_counts_control_reads_of_both_strands_together() {
    let snv = |read_id: usize, strand: Strand| -> Arc<VariantRecord> {
        let base: &str = if strand == Strand::Forward { "A" } else { "T" };
        Arc::new(VariantRecord::new(
            read_id,
            100,
            100,
            GraphOperation::new(
                0,
                7_674_224,
                strand.clone(),
                GraphOperationType::Downstream,
                0,
                7_674_226,
                strand,
                GraphOperationType::Upstream,
                base.into(),
                VariantType::SingleNucleotideVariant
            )
        ))
    };
    let a: Vec<Arc<VariantRecord>> = vec![snv(1, Strand::Forward)];
    let b: Vec<Arc<VariantRecord>> = vec![snv(11, Strand::Forward), snv(12, Strand::Reverse)];

    let a_only: Vec<Arc<VariantRecord>> = diff_dna_variant_records(a, b, 100_000, 1, 0.5, 0.5, 1_000, 0.01, 2, true, false);

    assert!(a_only.is_empty());
}


/// The control holds another base at the tumour SNV's position. Under the infinite sites
/// assumption a position mutates once, so the control SNV is taken to be the tumour's; without
/// it the tumour's base is another allele and stays. The same holds for an MNV.
#[test]
fn diff_dna_variant_records_applies_the_infinite_sites_assumption_to_substitutions() {
    let substitution = |read_id: usize, bases: &str| -> Arc<VariantRecord> {
        let variant_type: VariantType = if bases.len() == 1 {
            VariantType::SingleNucleotideVariant
        } else {
            VariantType::MultiNucleotideVariant
        };
        Arc::new(VariantRecord::new(
            read_id,
            100,
            100,
            GraphOperation::new(
                0,
                7_674_224,
                Strand::Forward,
                GraphOperationType::Downstream,
                0,
                7_674_225 + bases.len() as u32,
                Strand::Forward,
                GraphOperationType::Upstream,
                bases.into(),
                variant_type
            )
        ))
    };

    let kept = |a: &str, b: &str, apply_infinite_sites_assumption: bool| -> usize {
        diff_dna_variant_records(
            vec![substitution(1, a)],
            vec![substitution(11, b)],
            100_000, 1, 0.5, 0.5, 1_000, 0.01, 0,
            apply_infinite_sites_assumption,
            false
        ).len()
    };

    assert_eq!(kept("A", "T", true), 0);
    assert_eq!(kept("A", "T", false), 1);
    assert_eq!(kept("A", "A", false), 0);
    assert_eq!(kept("AC", "GT", true), 0);
    assert_eq!(kept("AC", "GT", false), 1);
    assert_eq!(kept("AC", "AC", false), 0);
}


/// A breakpoint held as split reads in one sample and as clips in the other is one variant: the
/// clip's breakend is one of the junction's, by chromosome and operation, and a clip has no size to
/// compare. The breakends may lie as far apart as an insertion of the clipped bases may (13 bases
/// for a 40-base clip, 3 for a 9-base one, at sequencing error 0.01). A clip that points the other
/// way, or lies farther away, is another variant: a normal read with a short clip 415 bases from a
/// breakend does not hold it. The infinite sites assumption is not needed.
#[test]
fn diff_dna_variant_records_matches_a_clip_with_a_split_read_at_either_breakend() {
    let split = |read_id: usize, chromosome_2: u16, position_2: u32| -> Arc<VariantRecord> {
        let variant_type: VariantType = if chromosome_2 == 0 { VariantType::Breakpoint } else { VariantType::Translocation };
        Arc::new(VariantRecord::new(
            read_id,
            100,
            101,
            GraphOperation::new(
                0,
                1_000,
                Strand::Forward,
                GraphOperationType::Downstream,
                chromosome_2,
                position_2,
                Strand::Forward,
                GraphOperationType::Upstream,
                "".into(),
                variant_type
            )
        ))
    };
    let clip = |read_id: usize, chromosome: u16, position: u32, operation: GraphOperationType, sequence: &str| -> Arc<VariantRecord> {
        Arc::new(VariantRecord::new(
            read_id,
            100,
            100 + sequence.len() as u32 - 1,
            GraphOperation::new(
                chromosome,
                position,
                Strand::Forward,
                operation,
                chromosome,
                position,
                Strand::Forward,
                GraphOperationType::Noop,
                sequence.into(),
                VariantType::Breakpoint
            )
        ))
    };
    let kept = |a: Arc<VariantRecord>, b: Arc<VariantRecord>, apply_infinite_sites_assumption: bool| -> usize {
        diff_dna_variant_records(vec![a], vec![b], 10_000, 1, 0.5, 0.5, 1_000, 0.01, 1, apply_infinite_sites_assumption, false).len()
    };
    let tail: &str = "ACGTACGTACGTACGTACGTAAAAACCCCCGGGGGTTTTT";
    let short_tail: &str = "ACGTAAACC";

    // The case holds the junction 1000 -> 5000 as a split read, the control as clips.
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_002, GraphOperationType::Downstream, tail), true), 0);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 5_012, GraphOperationType::Upstream, tail), true), 0);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_000, GraphOperationType::Downstream, tail), false), 0);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_002, GraphOperationType::Downstream, short_tail), true), 0);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_002, GraphOperationType::Upstream, tail), true), 1);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_014, GraphOperationType::Downstream, tail), true), 1);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 1_004, GraphOperationType::Downstream, short_tail), true), 1);
    assert_eq!(kept(split(1, 0, 5_000), clip(11, 0, 5_415, GraphOperationType::Upstream, short_tail), true), 1);

    // The case holds it as a clip, the control as a split read.
    assert_eq!(kept(clip(1, 0, 1_002, GraphOperationType::Downstream, tail), split(11, 0, 5_000), true), 0);
    assert_eq!(kept(clip(1, 0, 4_998, GraphOperationType::Upstream, tail), split(11, 0, 5_000), true), 0);

    // A translocation chr0:1000 -> chr1:5000, clipped in the control at either breakend.
    assert_eq!(kept(split(1, 1, 5_000), clip(11, 1, 5_002, GraphOperationType::Upstream, tail), true), 0);
    assert_eq!(kept(split(1, 1, 5_000), clip(11, 0, 998, GraphOperationType::Downstream, tail), true), 0);
    assert_eq!(kept(clip(1, 1, 5_002, GraphOperationType::Upstream, tail), split(11, 1, 5_000), true), 0);

    // A translocation and an intrachromosomal breakpoint that share a breakend stay two variants.
    assert_eq!(kept(split(1, 1, 5_000), split(11, 0, 3_000), true), 1);
}
