use crate::prelude::SpliceJunction;

use super::*;


/// The depth of a retained intron is its own reads: the reads retaining it and the reads
/// splicing it.
///
/// Each case is one intron of chr18, 1,000,101-1,000,900. One read retains it and the rest of its
/// reads splice it. A second intron, 2,000,101-2,000,900, is spliced by 1,000 reads and retained
/// by none: it takes no part in the depth of the first. A read that splices an intron retains
/// none, and holds no splicing event.
#[test]
fn rna_junction_read_counts_returns_retaining_reads_and_depth() {
    for num_reads in [2usize, 4, 300, 500] {
        let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
        // (intron start, reads retaining it, reads splicing it)
        for (intron_start, num_retaining_reads, num_splicing_reads) in [
            (1_000_101u32, 1usize, num_reads - 1),
            (2_000_101u32, 0usize, 1_000usize)
        ] {
            for _ in 0..num_retaining_reads {
                summaries.push(RNAReadCharacterizationSummary {
                    read_id: summaries.len(),
                    chromosome: Some(1),
                    reference_start: intron_start - 100,
                    reference_end: intron_start + 899,
                    ends_at_reference_start: true,
                    ends_at_reference_end: true,
                    splice_junctions: Vec::new(),
                    exons: vec![(1, intron_start - 100, intron_start + 899)],
                    offset: 0,
                    length: 0
                });
            }
            for _ in 0..num_splicing_reads {
                summaries.push(RNAReadCharacterizationSummary {
                    read_id: summaries.len(),
                    chromosome: Some(1),
                    reference_start: intron_start - 100,
                    reference_end: intron_start + 899,
                    ends_at_reference_start: true,
                    ends_at_reference_end: true,
                    splice_junctions: vec![
                        SpliceJunction::new(1, 1, intron_start, intron_start + 799, Strand::Forward, Strand::Forward)
                    ],
                    exons: vec![(1, intron_start - 100, intron_start - 1), (1, intron_start + 800, intron_start + 899)],
                    offset: 0,
                    length: 0
                });
            }
        }

        let read_counts: RNAJunctionReadCounts = RNAJunctionReadCounts::new(&summaries, 2);

        assert_eq!(read_counts.get_read_counts(&summaries[0]), vec![(1, num_reads as ReadDepth)], "{} reads", num_reads);
        for summary in summaries[1..].iter() {
            assert!(read_counts.get_read_counts(summary).is_empty(), "{} reads, read {}", num_reads, summary.read_id);
        }
        assert_eq!(read_counts.get_depths(), HashSet::from([num_reads as ReadDepth, 1_000]), "{} reads", num_reads);
    }
}


/// A read that holds the exonic base on one side of an intron and ends among the bases of the
/// intron retains it. Where its block stops there and the read goes on, it says nothing.
///
/// The reads sit on three exons of TP53 (chr17, reverse strand) and the two introns between
/// them, 7,674,972-7,675,052 and 7,675,237-7,675,993:
///
///   302 reads  splice both introns
///    40 reads  retain the first intron and splice the second
///     8 reads  splice the first intron, run 364 bases into the second and stop
///
/// Where the 8 reads end there, they retain the second intron: 8 retaining reads of a depth of
/// 350, with the 342 that splice it. The first intron is retained by the 40, of the 350 with the
/// 310 that splice it. Where the 8 reads go on, no read retains the second intron, and its depth
/// is the 342 that splice it.
#[test]
fn rna_junction_read_counts_returns_reads_ending_inside_intron() {
    // (the 8 reads end at the last base of their blocks, counts of each of the 8, depths)
    let cases: Vec<(bool, Vec<(ReadSupport, ReadDepth)>, HashSet<ReadDepth>)> = vec![
        (true, vec![(8, 350)], HashSet::from([350])),
        (false, Vec::new(), HashSet::from([350, 342]))
    ];
    for (ends_at_reference_end, expected_read_counts, expected_depths) in cases {
        let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
        for read_id in 0..302usize {
            summaries.push(RNAReadCharacterizationSummary {
                read_id,
                chromosome: Some(0),
                reference_start: 7_674_859,
                reference_end: 7_676_272,
                ends_at_reference_start: true,
                ends_at_reference_end: true,
                splice_junctions: vec![
                    SpliceJunction::new(0, 0, 7_675_993, 7_675_237, Strand::Reverse, Strand::Reverse),
                    SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
                ],
                exons: vec![(0, 7_675_994, 7_676_272), (0, 7_675_053, 7_675_236), (0, 7_674_859, 7_674_971)],
                offset: 0,
                length: 0
            });
        }
        for read_id in 302..342usize {
            summaries.push(RNAReadCharacterizationSummary {
                read_id,
                chromosome: Some(0),
                reference_start: 7_674_859,
                reference_end: 7_676_272,
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
        for read_id in 342..350usize {
            summaries.push(RNAReadCharacterizationSummary {
                read_id,
                chromosome: Some(0),
                reference_start: 7_674_859,
                reference_end: 7_675_600,
                ends_at_reference_start: true,
                ends_at_reference_end,
                splice_junctions: vec![
                    SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
                ],
                exons: vec![(0, 7_675_053, 7_675_600), (0, 7_674_859, 7_674_971)],
                offset: 0,
                length: 0
            });
        }

        let read_counts: RNAJunctionReadCounts = RNAJunctionReadCounts::new(&summaries, 2);

        for summary in summaries.iter() {
            let expected: Vec<(ReadSupport, ReadDepth)> = match summary.read_id {
                0..=301 => Vec::new(),
                302..=341 => vec![(40, 350)],
                _ => expected_read_counts.clone()
            };
            assert_eq!(
                read_counts.get_read_counts(summary),
                expected,
                "read {}, the 8 reads end at the last base of their blocks: {}",
                summary.read_id, ends_at_reference_end
            );
        }
        assert_eq!(
            read_counts.get_depths(),
            expected_depths,
            "the 8 reads end at the last base of their blocks: {}",
            ends_at_reference_end
        );
    }
}
