use super::*;


/// The depth of a novel splice junction is the reads splicing it or an intron overlapping it.
///
/// Each case holds an annotated intron of chr18, 1,000,101-1,000,900, and a novel splice
/// junction overlapping it, 1,000,105-1,000,900. One read splices the novel splice junction and
/// the rest of its reads splice the annotated intron, so its depth is the two together. A second
/// annotated intron, 2,000,101-2,000,900, is spliced by 1,000 reads: it overlaps nothing and
/// takes no part in the depth. A read that splices only annotated introns holds no splicing
/// event.
#[test]
fn rna_novel_junction_read_counts_returns_splicing_reads_and_depth() {
    let annotated_introns: HashSet<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = HashSet::from([
        (1, 1_000_101, 1_000_900),
        (1, 2_000_101, 2_000_900)
    ]);

    for num_reads in [2usize, 4, 300, 500] {
        let mut splice_junctions: Vec<Vec<SpliceJunction>> = Vec::new();
        // (first intron base, number of reads)
        for (intron_start, num_splicing_reads) in [
            (1_000_105u32, 1usize),
            (1_000_101u32, num_reads - 1),
            (2_000_101u32, 1_000usize)
        ] {
            for _ in 0..num_splicing_reads {
                splice_junctions.push(vec![
                    SpliceJunction::new(1, 1, intron_start, intron_start / 1_000_000 * 1_000_000 + 900, Strand::Forward, Strand::Forward)
                ]);
            }
        }

        let read_counts: RNANovelJunctionReadCounts = RNANovelJunctionReadCounts::new(
            &splice_junctions,
            |intron| annotated_introns.contains(intron)
        );

        assert_eq!(read_counts.get_read_counts(&splice_junctions[0]), vec![(1, num_reads as ReadDepth)], "{} reads", num_reads);
        for read_splice_junctions in splice_junctions[1..].iter() {
            assert!(read_counts.get_read_counts(read_splice_junctions).is_empty(), "{} reads", num_reads);
        }
        assert_eq!(read_counts.get_depths(), HashSet::from([num_reads as ReadDepth]), "{} reads", num_reads);
    }
}


/// A read counts once at a novel splice junction, however many of the introns overlapping it it
/// splices and however many times its chain spells one. A junction joining two chromosomes is
/// not an intron.
///
/// On chr17 the annotated intron 1,000,101-1,000,900 and two novel splice junctions,
/// 1,000,105-1,000,900 and 1,000,101-1,000,950, all overlap one another:
///
///    5 reads  splice 1,000,105-1,000,900 on two laps of a circle
///    3 reads  splice 1,000,101-1,000,950 and the annotated intron
///   10 reads  splice the annotated intron
///    2 reads  join chr17 1,000,105 to chr18 1,000,900
///
/// Each novel splice junction has a depth of the 18 reads on chr17.
#[test]
fn rna_novel_junction_read_counts_counts_a_read_once() {
    let mut splice_junctions: Vec<Vec<SpliceJunction>> = Vec::new();
    for _ in 0..5 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 0, 1_000_105, 1_000_900, Strand::Forward, Strand::Forward),
            SpliceJunction::new(0, 0, 1_000_105, 1_000_900, Strand::Forward, Strand::Forward)
        ]);
    }
    for _ in 0..3 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 0, 1_000_101, 1_000_950, Strand::Forward, Strand::Forward),
            SpliceJunction::new(0, 0, 1_000_101, 1_000_900, Strand::Forward, Strand::Forward)
        ]);
    }
    for _ in 0..10 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 0, 1_000_101, 1_000_900, Strand::Forward, Strand::Forward)
        ]);
    }
    for _ in 0..2 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 1, 1_000_105, 1_000_900, Strand::Forward, Strand::Forward)
        ]);
    }

    let read_counts: RNANovelJunctionReadCounts = RNANovelJunctionReadCounts::new(
        &splice_junctions,
        |intron| *intron == (0, 1_000_101, 1_000_900)
    );

    assert_eq!(read_counts.get_read_counts(&splice_junctions[0]), vec![(5, 18)]);
    assert_eq!(read_counts.get_read_counts(&splice_junctions[5]), vec![(3, 18)]);
    assert!(read_counts.get_read_counts(&splice_junctions[8]).is_empty());
    assert!(read_counts.get_read_counts(&splice_junctions[18]).is_empty());
    assert_eq!(read_counts.get_depths(), HashSet::from([18]));
}
