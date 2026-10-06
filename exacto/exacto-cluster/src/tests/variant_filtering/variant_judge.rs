use std::collections::HashSet;

use super::*;


/// A call needs the larger floor of its two sides, each read at its own depth, raised in a
/// repeat. A side no read of the cluster covers has no floor a call could clear.
///
/// At a sequencing error and slippage rate of 0.01 and a max FPR of 1e-6:
///
///   Repeat length   Depths      Minimum reads
///   0               30, 30      9
///   0               2, 100      19, the deeper side's
///   0               100, 2      19
///   10              30, 30      27
///   0               0, 30       none (ReadSupport::MAX)
///   0               30, 0       none
#[test]
fn min_reads_at_site_returns_larger_floor_of_two_sides() {
    let depths: HashSet<u32> = [2u32, 30, 100].into_iter().collect();
    let read_support_index: RNAVariantReadSupportIndex = RNAVariantReadSupportIndex::new(
        &depths,
        30,
        0.01,
        0.01,
        1e-6,
        1
    );

    assert_eq!(min_reads_at_site(&read_support_index, 0, (30, 30)), 9);
    assert_eq!(min_reads_at_site(&read_support_index, 0, (2, 100)), 19);
    assert_eq!(min_reads_at_site(&read_support_index, 0, (100, 2)), 19);
    assert_eq!(min_reads_at_site(&read_support_index, 10, (30, 30)), 27);
    assert_eq!(min_reads_at_site(&read_support_index, 0, (0, 30)), ReadSupport::MAX);
    assert_eq!(min_reads_at_site(&read_support_index, 0, (30, 0)), ReadSupport::MAX);
}
