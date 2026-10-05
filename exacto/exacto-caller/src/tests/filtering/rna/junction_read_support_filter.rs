use std::collections::HashSet;

use super::*;


/// A retained intron is held to the shared splicing event error model at its own depth, and the
/// reads retaining one that falls short fail.
///
/// The reads sit on three exons of TP53 (chr17, reverse strand) and the two introns between
/// them, 7,674,972-7,675,052 and 7,675,237-7,675,993:
///
///    14 reads  retain the first intron and splice the second
///   286 reads  splice the first intron and retain the second
///
/// Each intron is retained or spliced by all 300 reads, a minimum of 46 at an error rate of 0.01
/// and a max FPR of 1e-6. The first is retained by 14 reads, which fail. The second is retained
/// by 286, which pass.
#[test]
fn rna_junction_read_support_filter_fails_reads_retaining_intron_below_min_read_support() {
    let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
    for read_id in 0..14usize {
        summaries.push(RNAReadCharacterizationSummary {
            read_id,
            chromosome: Some(0),
            reference_start: 7_676_272,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_993, 7_675_237, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_994, 7_676_272), (0, 7_674_859, 7_675_236)],
            offset: 0,
            length: 0
        });
    }
    for read_id in 14..300usize {
        summaries.push(RNAReadCharacterizationSummary {
            read_id,
            chromosome: Some(0),
            reference_start: 7_676_272,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_053, 7_676_272), (0, 7_674_859, 7_674_971)],
            offset: 0,
            length: 0
        });
    }

    let read_counts: RNAJunctionReadCounts = RNAJunctionReadCounts::new(&summaries, 2);
    let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &read_counts.get_depths(),
        0.01,
        1e-6,
        2
    );
    let read_filter: RNAJunctionReadSupportFilter<RNAJunctionReadCounts> = RNAJunctionReadSupportFilter::new(
        &read_counts,
        &read_support_index
    );
    let failed_read_ids: HashSet<usize> = summaries
        .iter()
        .filter(|summary| !read_filter.passes(summary))
        .map(|summary| summary.read_id)
        .collect();

    assert_eq!(read_counts.get_depths(), HashSet::from([300]));
    assert_eq!(read_support_index.get_min_read_support(300), 46);
    assert_eq!(failed_read_ids, (0..14).collect::<HashSet<usize>>());
}


/// A novel splice junction is held to the shared splicing event error model at its own depth,
/// and the reads splicing one that falls short fail.
///
/// The reads sit on two introns of TP53 (chr17, reverse strand), 7,676,404-7,676,520 and
/// 7,676,623-7,687,376. All 300 splice the first. At the second, which no read splices as the
/// annotation holds it:
///
///   286 reads  splice 7,676,623-7,687,386, the second intron and the 10 bases after it
///    14 reads  splice 7,676,623-7,687,377, the second intron and the base after it
///
/// Each novel splice junction overlaps the other, so each has a depth of 300 reads and a minimum
/// of 46 at an error rate of 0.01 and a max FPR of 1e-6. The 14 reads fail and the 286 pass.
#[test]
fn rna_junction_read_support_filter_fails_reads_splicing_novel_junction_below_min_read_support() {
    let annotated_introns: HashSet<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = HashSet::from([
        (0, 7_676_404, 7_676_520),
        (0, 7_676_623, 7_687_376)
    ]);

    let mut splice_junctions: Vec<Vec<SpliceJunction>> = Vec::new();
    for _ in 0..286 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 0, 7_687_386, 7_676_623, Strand::Reverse, Strand::Reverse),
            SpliceJunction::new(0, 0, 7_676_520, 7_676_404, Strand::Reverse, Strand::Reverse)
        ]);
    }
    for _ in 286..300 {
        splice_junctions.push(vec![
            SpliceJunction::new(0, 0, 7_687_377, 7_676_623, Strand::Reverse, Strand::Reverse),
            SpliceJunction::new(0, 0, 7_676_520, 7_676_404, Strand::Reverse, Strand::Reverse)
        ]);
    }

    let read_counts: RNANovelJunctionReadCounts = RNANovelJunctionReadCounts::new(
        &splice_junctions,
        |intron| annotated_introns.contains(intron)
    );
    let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &read_counts.get_depths(),
        0.01,
        1e-6,
        2
    );
    let read_filter: RNAJunctionReadSupportFilter<RNANovelJunctionReadCounts> = RNAJunctionReadSupportFilter::new(
        &read_counts,
        &read_support_index
    );
    let failed_reads: HashSet<usize> = splice_junctions
        .iter()
        .enumerate()
        .filter(|(_, read_splice_junctions)| !read_filter.passes(read_splice_junctions))
        .map(|(index, _)| index)
        .collect();

    assert_eq!(read_counts.get_depths(), HashSet::from([300]));
    assert_eq!(read_support_index.get_min_read_support(300), 46);
    assert_eq!(failed_reads, (286..300).collect::<HashSet<usize>>());
}


/// A read passes at exactly the minimum read support of its event, and fails one read short.
///
/// The reads sit on the two TP53 introns of the retention test above. k reads retain the first
/// intron and splice the second, and the other 300 - k splice the first and retain the second,
/// so each intron has a depth of 300 and a minimum of 46.
#[test]
fn rna_junction_read_support_filter_passes_reads_at_min_read_support() {
    // (reads retaining the first intron, reads that fail)
    for (num_retaining_reads, expected_failed_read_ids) in [(46usize, HashSet::new()), (45, (0..45).collect::<HashSet<usize>>())] {
        let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
        for read_id in 0..num_retaining_reads {
            summaries.push(RNAReadCharacterizationSummary {
                read_id,
                chromosome: Some(0),
                reference_start: 7_676_272,
                reference_end: 7_674_859,
                ends_at_reference_start: true,
                ends_at_reference_end: true,
                splice_junctions: vec![
                    SpliceJunction::new(0, 0, 7_675_993, 7_675_237, Strand::Reverse, Strand::Reverse)
                ],
                exons: vec![(0, 7_675_994, 7_676_272), (0, 7_674_859, 7_675_236)],
                offset: 0,
                length: 0
            });
        }
        for read_id in num_retaining_reads..300 {
            summaries.push(RNAReadCharacterizationSummary {
                read_id,
                chromosome: Some(0),
                reference_start: 7_676_272,
                reference_end: 7_674_859,
                ends_at_reference_start: true,
                ends_at_reference_end: true,
                splice_junctions: vec![
                    SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
                ],
                exons: vec![(0, 7_675_053, 7_676_272), (0, 7_674_859, 7_674_971)],
                offset: 0,
                length: 0
            });
        }

        let read_counts: RNAJunctionReadCounts = RNAJunctionReadCounts::new(&summaries, 2);
        let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
            &read_counts.get_depths(),
            0.01,
            1e-6,
            2
        );
        let read_filter: RNAJunctionReadSupportFilter<RNAJunctionReadCounts> = RNAJunctionReadSupportFilter::new(
            &read_counts,
            &read_support_index
        );
        let failed_read_ids: HashSet<usize> = summaries
            .iter()
            .filter(|summary| !read_filter.passes(summary))
            .map(|summary| summary.read_id)
            .collect();

        assert_eq!(read_support_index.get_min_read_support(300), 46);
        assert_eq!(failed_read_ids, expected_failed_read_ids, "{} reads retain the first intron", num_retaining_reads);
    }
}
