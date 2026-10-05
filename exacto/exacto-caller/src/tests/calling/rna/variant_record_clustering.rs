use exacto_core::prelude::*;

use super::*;


#[test]
fn spliced_distance_charges_only_exonic_bases_outside_and_across_the_span() {
    // Introns of one chain: exon A | 100..=199 | exon B | 300..=399 | exon C
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 100, 199, Strand::Forward, Strand::Forward),
        SpliceJunction::new(0, 0, 300, 399, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);

    assert_eq!(ladders.spliced_distance(0, 10, 50), 40);      // both before the first intron
    assert_eq!(ladders.spliced_distance(0, 450, 900), 450);   // both beyond the last intron
    assert_eq!(ladders.spliced_distance(0, 250, 450), 100);   // 250..=299 in B (50) + 400..=450 in C (50)
    assert_eq!(ladders.spliced_distance(0, 90, 210), 20);     // 91..=99 in A (9) + 200..=210 in B (11)
    assert_eq!(ladders.spliced_distance(1, 10, 50), 40);      // chromosome with no ladder: genomic
}
/// One junction spelled at jittered breakends must cluster as ONE junction, and a junction
/// whose side 1 lies beyond `max_breakend_distance` must keep its own. Each side is pooled
/// on its own, so the distant junction shares the side-2 pool with the others and still
/// separates on side 1. The distance is spliced: with an intron between the near and the
/// distant side-1 breakends, all three spellings pool after all.
#[test]
fn pool_breakpoints_pools_jittered_breakends_and_keeps_distant_ones_apart() {
    let breakpoint = |position_1: u32, position_2: u32| -> GraphOperation {
        GraphOperation::new(
            0,
            position_1,
            Strand::Forward,
            GraphOperationType::Downstream,
            1,
            position_2,
            Strand::Forward,
            GraphOperationType::Upstream,
            "".into(),
            VariantType::Breakpoint
        )
    };

    // Records in read id order, which is the order the clusterer sorts them into.
    let major: GraphOperation = breakpoint(7_000_000, 3_000_000);   // reads 1..=120
    let minor: GraphOperation = breakpoint(7_000_004, 3_000_003);   // reads 121..=200, within tolerance on both sides
    let distant: GraphOperation = breakpoint(7_000_200, 3_000_000); // reads 201..=230, side 1 beyond it
    let variant_records: Vec<Arc<VariantRecord>> = (1..=230)
        .map(|read_id| {
            let spelling: GraphOperation = if read_id <= 120 {
                major.clone()
            } else if read_id <= 200 {
                minor.clone()
            } else {
                distant.clone()
            };
            Arc::new(VariantRecord::new(read_id, 100, 100, spelling))
        })
        .collect();

    // No introns: distances are genomic. Side 1 splits the distant junction off, side 2
    // holds all three; pool ids follow each pool's smallest member.
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&[]);
    let pool_1: Vec<usize> = pool_breakpoints(&variant_records, Side::One, &ladders, 50);
    let pool_2: Vec<usize> = pool_breakpoints(&variant_records, Side::Two, &ladders, 50);
    assert!(pool_1[..200].iter().all(|&pool| pool == 0));
    assert!(pool_1[200..].iter().all(|&pool| pool == 1));
    assert!(pool_2.iter().all(|&pool| pool == 0));

    // A junction is a pair of pools: two clusters, the jittered spellings together.
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_rna_variant_records(variant_records.clone(), &ladders, 50);
    assert_eq!(clusters.len(), 2);
    let read_ids = |cluster: &VariantRecordCluster| -> HashSet<usize> {
        cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect()
    };
    assert_eq!(read_ids(&clusters[0]), (1..=200).collect::<HashSet<usize>>());
    assert_eq!(read_ids(&clusters[1]), (201..=230).collect::<HashSet<usize>>());

    // An intron between the near and the distant side-1 breakends: 7,000,010..=7,000,190 is
    // skipped, so 7,000,000 and 7,000,200 sit 19 spliced bp apart and everything pools.
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 7_000_010, 7_000_190, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);
    assert_eq!(ladders.spliced_distance(0, 7_000_000, 7_000_200), 19);
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_rna_variant_records(variant_records, &ladders, 50);
    assert_eq!(clusters.len(), 1);
    assert_eq!(read_ids(&clusters[0]), (1..=230).collect::<HashSet<usize>>());
}


/// One deletion event spelled two ways must cluster as ONE cluster, the H1-2 case: aligner
/// jitter in a repeat spells the same excision as 15 bp in some reads and 16 bp in others.
/// The size gate still separates a spelling that is not the same event: a 40-bp deletion
/// from the same left anchor is under half the size proportion and keeps its own cluster.
/// The call built from the pooled cluster then carries every carrier and the modal spelling,
/// so neither spelling's carriers count as reference evidence against the other.
#[test]
fn cluster_non_breakpoint_rna_variant_records_pools_deletion_spellings() {
    // Flanking-anchor grammar: the deleted bases sit on `[position_1 + 1, position_2 - 1]`,
    // so a 15-bp deletion starting at `p + 1` anchors at `p` and `p + 16`.
    let deletion = |position_1: u32, position_2: u32| -> GraphOperation {
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
    };

    let major: GraphOperation = deletion(6_000_100, 6_000_116); // 15-bp spelling, reads 1..=120
    let minor: GraphOperation = deletion(6_000_100, 6_000_117); // 16-bp spelling, reads 121..=200
    let other: GraphOperation = deletion(6_000_100, 6_000_141); // 40-bp spelling, reads 201..=230
    let variant_records: Vec<Arc<VariantRecord>> = (1..=230)
        .map(|read_id| {
            let spelling: GraphOperation = if read_id <= 120 {
                major.clone()
            } else if read_id <= 200 {
                minor.clone()
            } else {
                other.clone()
            };
            Arc::new(VariantRecord::new(read_id, 100, 100, spelling))
        })
        .collect();

    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&[]);
    let clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_rna_variant_records(
        &variant_records,
        &ladders,
        0.5,        // min_size_proportion
        0.5,        // max_ins_norm_edit_distance
        50,         // max_clustering_distance
        0.01        // sequencing_error
    );

    // Two clusters: the two jittered spellings together, the 40-bp deletion alone.
    assert_eq!(clusters.len(), 2);
    let read_ids = |cluster: &VariantRecordCluster| -> HashSet<usize> {
        cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect()
    };
    assert_eq!(read_ids(&clusters[0]), (1..=200).collect::<HashSet<usize>>());
    assert_eq!(read_ids(&clusters[1]), (201..=230).collect::<HashSet<usize>>());

    // The call made from the pooled cluster: every carrier is a member, at the modal spelling.
    let variant_call: VariantCall = VariantCall::from_variant_records(
        0,
        clusters[0].get_variant_records().iter().map(Arc::as_ref).cloned().collect(),
        0,
        4,
        6,
        2
    );
    assert_eq!(variant_call.get_read_ids().len(), 200);
    assert_eq!(*variant_call.get_consensus_graph_operation(), major);
}

/// A position inside an intron is a position on a read that retains the intron, so the bases
/// of the intron between two positions count unless the intron lies whole between them.
///
/// Introns of one chain: exon A | 100..=199 | exon B | 300..=399 | exon C
///
///   Positions   Lie                                  Introns whole between   Distance
///   120, 180    inside the first intron              none                    60
///   90, 150     on exon A, inside the first intron   none                    60
///   150, 250    inside the first intron, on exon B   none                    100
///   150, 350    inside the first, inside the second  none                    200
///   150, 450    inside the first intron, on exon C   the second              200
///   99, 200     on exon A, on exon B                 the first               1
///   90, 450     on exon A, on exon C                 both                    160
#[test]
fn spliced_distance_returns_genomic_distance_inside_intron() {
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 100, 199, Strand::Forward, Strand::Forward),
        SpliceJunction::new(0, 0, 300, 399, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);

    // (position, position, distance)
    let cases: Vec<(u32, u32, u32)> = vec![
        (120, 180, 60),
        (90, 150, 60),
        (150, 250, 100),
        (150, 350, 200),
        (150, 450, 200),
        (99, 200, 1),
        (90, 450, 160)
    ];
    for (a, b, expected) in cases {
        assert_eq!(ladders.spliced_distance(0, a, b), expected, "{} to {}", a, b);
        assert_eq!(ladders.spliced_distance(0, b, a), expected, "{} to {}", b, a);
    }
}


/// A sweep over positions in genomic order goes on past an intron when a position inside it is
/// out of tolerance, and is over when a position on an exon is.
///
/// Introns of one chain: exon A | 100..=199 | exon B | 300..=399 | exon C
/// Positions, in order: 90, 120, 150, 180, 205, 260, 350, 420
///
///   From   Position out of tolerance   Lies                       Sweep goes on from
///   90     150                         inside the first intron    205
///   120    180                         inside the intron of 120   nowhere
///   90     260                         on exon B                  nowhere
///   205    350                         inside the second intron   420
///   180    350                         inside the second intron   420
#[test]
fn next_sweep_index_returns_first_position_past_intron() {
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 100, 199, Strand::Forward, Strand::Forward),
        SpliceJunction::new(0, 0, 300, 399, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);
    let positions: Vec<u32> = vec![90, 120, 150, 180, 205, 260, 350, 420];

    // (from, index of the position out of tolerance, index the sweep goes on from)
    let cases: Vec<(u32, usize, Option<usize>)> = vec![
        (90, 2, Some(4)),
        (120, 3, None),
        (90, 5, None),
        (205, 6, Some(7)),
        (180, 6, Some(7))
    ];
    for (a, index, expected) in cases {
        assert_eq!(
            ladders.next_sweep_index(0, a, &positions, index),
            expected,
            "from {}, position {}",
            a,
            positions[index]
        );
    }

    // A chromosome with no ladder has no intron to go past.
    assert_eq!(ladders.next_sweep_index(1, 90, &positions, 2), None);
}


/// A breakend inside a retained intron pools with nothing across the intron, and does not keep
/// the breakends on the two exons beside the intron from pooling.
///
/// The intron 7,000,010..=7,000,190 lies between the breakends of reads 1..=40 and 51..=90, so
/// they sit 19 spliced bp apart. Reads 41..=50 retain the intron and break inside it, 100 bp
/// from both.
#[test]
fn pool_breakpoints_keeps_breakends_inside_intron_apart() {
    let breakpoint = |position_1: u32| -> GraphOperation {
        GraphOperation::new(
            0,
            position_1,
            Strand::Forward,
            GraphOperationType::Downstream,
            1,
            3_000_000,
            Strand::Forward,
            GraphOperationType::Upstream,
            "".into(),
            VariantType::Breakpoint
        )
    };
    let variant_records: Vec<Arc<VariantRecord>> = (1..=90)
        .map(|read_id| {
            let spelling: GraphOperation = if read_id <= 40 {
                breakpoint(7_000_000)
            } else if read_id <= 50 {
                breakpoint(7_000_100)
            } else {
                breakpoint(7_000_200)
            };
            Arc::new(VariantRecord::new(read_id, 100, 100, spelling))
        })
        .collect();
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 7_000_010, 7_000_190, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);

    let pool_1: Vec<usize> = pool_breakpoints(&variant_records, Side::One, &ladders, 50);
    assert!(pool_1[..40].iter().all(|&pool| pool == 0));
    assert!(pool_1[40..50].iter().all(|&pool| pool == 1));
    assert!(pool_1[50..].iter().all(|&pool| pool == 0));
}


/// The records of reads retaining an intron pool by where they sit in the intron.
///
/// The intron 6,070,653..=6,080,370 is an intron of the chain. Reads 1..=9 retain it. Reads
/// 1..=3 carry one deletion of 1 base at 6,077,089; every read carries a deletion of 1 base of
/// its own as well, 1,000 bases apart from the next.
#[test]
fn cluster_non_breakpoint_rna_variant_records_keeps_records_inside_intron_apart() {
    let deletion = |position: u32| -> GraphOperation {
        GraphOperation::new(
            0,
            position - 1,
            Strand::Forward,
            GraphOperationType::Downstream,
            0,
            position + 1,
            Strand::Forward,
            GraphOperationType::Upstream,
            "".into(),
            VariantType::Deletion
        )
    };
    let mut variant_records: Vec<Arc<VariantRecord>> = Vec::new();
    for read_id in 1..=9usize {
        if read_id <= 3 {
            variant_records.push(Arc::new(VariantRecord::new(read_id, 6_437, 6_438, deletion(6_077_089))));
        }
        variant_records.push(Arc::new(VariantRecord::new(read_id, 500, 501, deletion(6_070_700 + 1_000 * read_id as u32))));
    }
    let splice_junctions: Vec<SpliceJunction> = vec![
        SpliceJunction::new(0, 0, 6_070_653, 6_080_370, Strand::Forward, Strand::Forward)
    ];
    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&splice_junctions);

    let clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_rna_variant_records(
        &variant_records,
        &ladders,
        0.5,        // min_size_proportion
        0.5,        // max_ins_norm_edit_distance
        1_000,      // max_clustering_distance
        0.01        // sequencing_error
    );

    // One cluster of the three reads sharing a deletion, and one per deletion of a read's own.
    assert_eq!(clusters.len(), 10);
    let mut read_ids: Vec<Vec<usize>> = clusters
        .iter()
        .map(|cluster| {
            let mut read_ids: Vec<usize> = cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect();
            read_ids.sort_unstable();
            read_ids
        })
        .collect();
    read_ids.sort();
    assert_eq!(
        read_ids,
        vec![vec![1], vec![1, 2, 3], vec![2], vec![3], vec![4], vec![5], vec![6], vec![7], vec![8], vec![9]]
    );
}


/// Three reads hold A, two hold T and one holds C at one position: three clusters, one per base.
#[test]
fn cluster_non_breakpoint_rna_variant_records_keeps_snv_bases_at_one_position_apart() {
    let variant_records: Vec<Arc<VariantRecord>> = [(1, "A"), (2, "A"), (3, "A"), (4, "T"), (5, "T"), (6, "C")]
        .iter()
        .map(|&(read_id, base)| {
            Arc::new(VariantRecord::new(
                read_id,
                100,
                100,
                GraphOperation::new(
                    0,
                    6_000_100,
                    Strand::Forward,
                    GraphOperationType::Downstream,
                    0,
                    6_000_102,
                    Strand::Forward,
                    GraphOperationType::Upstream,
                    base.into(),
                    VariantType::SingleNucleotideVariant
                )
            ))
        })
        .collect();

    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&[]);
    let clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_rna_variant_records(
        &variant_records,
        &ladders,
        0.5,        // min_size_proportion
        0.5,        // max_ins_norm_edit_distance
        50,         // max_clustering_distance
        0.01        // sequencing_error
    );

    let mut alleles: Vec<(String, Vec<usize>)> = clusters
        .iter()
        .map(|cluster| {
            let records: Vec<&Arc<VariantRecord>> = cluster.get_variant_records().iter().collect();
            let mut read_ids: Vec<usize> = records.iter().map(|record| record.get_read_id()).collect();
            read_ids.sort();
            (records[0].get_standardized_sequence(), read_ids)
        })
        .collect();
    alleles.sort();
    assert_eq!(
        alleles,
        vec![("A".to_string(), vec![1, 2, 3]), ("C".to_string(), vec![6]), ("T".to_string(), vec![4, 5])]
    );
}


/// The reads of one transcript sequenced from either end spell a substitution on both strands:
/// three reads aligned forward hold A, three aligned reverse hold T in read orientation, which is
/// A on the forward strand. One cluster of six (exacto-cluster review F-06).
#[test]
fn cluster_non_breakpoint_rna_variant_records_pools_snv_across_alignment_strands() {
    let variant_records: Vec<Arc<VariantRecord>> = [(1, Strand::Forward, "A"), (2, Strand::Forward, "A"), (3, Strand::Forward, "A"),
        (4, Strand::Reverse, "T"), (5, Strand::Reverse, "T"), (6, Strand::Reverse, "T")]
        .iter()
        .map(|(read_id, strand, base)| {
            Arc::new(VariantRecord::new(
                *read_id,
                100,
                100,
                GraphOperation::new(
                    0,
                    6_000_100,
                    strand.clone(),
                    GraphOperationType::Downstream,
                    0,
                    6_000_102,
                    strand.clone(),
                    GraphOperationType::Upstream,
                    (*base).into(),
                    VariantType::SingleNucleotideVariant
                )
            ))
        })
        .collect();

    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&[]);
    let clusters: Vec<VariantRecordCluster> = cluster_non_breakpoint_rna_variant_records(
        &variant_records,
        &ladders,
        0.5,        // min_size_proportion
        0.5,        // max_ins_norm_edit_distance
        50,         // max_clustering_distance
        0.01        // sequencing_error
    );

    let mut read_ids: Vec<Vec<usize>> = clusters
        .iter()
        .map(|cluster| {
            let mut read_ids: Vec<usize> = cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect();
            read_ids.sort();
            read_ids
        })
        .collect();
    read_ids.sort();
    assert_eq!(read_ids, vec![vec![1, 2, 3, 4, 5, 6]]);
}


/// The reads of one fusion sequenced from either end spell its junction at the same two bases
/// with the same operation types, both strands flipped and the untemplated bases reverse
/// complemented. One cluster of six, not one per orientation: genotyping takes a call's members
/// as its carriers, so two calls would split the fusion's cluster by read orientation.
#[test]
fn cluster_breakpoint_rna_variant_records_pools_breakends_across_alignment_strands() {
    let variant_records: Vec<Arc<VariantRecord>> = [(1, Strand::Forward, "AGG"), (2, Strand::Forward, "AGG"), (3, Strand::Forward, "AGG"),
        (4, Strand::Reverse, "CCT"), (5, Strand::Reverse, "CCT"), (6, Strand::Reverse, "CCT")]
        .iter()
        .map(|(read_id, strand, sequence)| {
            Arc::new(VariantRecord::new(
                *read_id,
                100,
                104,
                GraphOperation::new(
                    0,
                    3_489_340,
                    strand.clone(),
                    GraphOperationType::Downstream,
                    0,
                    6_087_991,
                    strand.clone(),
                    GraphOperationType::Upstream,
                    (*sequence).into(),
                    VariantType::FusionGene
                )
            ))
        })
        .collect();

    let ladders: ExonLadders = ExonLadders::from_splice_junctions(&[]);
    let clusters: Vec<VariantRecordCluster> = cluster_breakpoint_rna_variant_records(variant_records, &ladders, 50);

    let mut read_ids: Vec<Vec<usize>> = clusters
        .iter()
        .map(|cluster| {
            let mut read_ids: Vec<usize> = cluster.get_variant_records().iter().map(|record| record.get_read_id()).collect();
            read_ids.sort();
            read_ids
        })
        .collect();
    read_ids.sort();
    assert_eq!(read_ids, vec![vec![1, 2, 3, 4, 5, 6]]);
}
