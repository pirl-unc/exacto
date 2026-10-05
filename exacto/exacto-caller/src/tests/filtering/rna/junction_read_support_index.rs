use rayon::prelude::*;
use rayon::ThreadPool;

use super::*;


/// Minimum read supports against the upper tail of BetaBinomial(n, 0.99, 98.01), an error rate
/// and an overdispersion of 0.01, summed term by term from SciPy's PMF (never as one minus a CDF
/// near one).
#[test]
fn rna_junction_read_support_index_returns_matches() {
    // (max FPR, depth, minimum read support)
    for (max_fpr, depth, min_read_support) in [
        (1e-6f64, 0u32, 2u32),
        (1e-6, 1, 2),
        // No count at this depth is significant: the minimum is out of reach.
        (1e-6, 2, 3),
        (1e-6, 4, 4),
        (1e-6, 10, 6),
        (1e-6, 20, 8),
        (1e-6, 50, 12),
        (1e-6, 100, 19),
        (1e-6, 145, 25),
        (1e-6, 200, 33),
        (1e-6, 292, 45),
        (1e-6, 300, 46),
        (1e-6, 500, 72),
        (1e-4, 300, 32)
    ] {
        let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
            &HashSet::from([depth]),
            0.01,
            max_fpr,
            2
        );

        assert_eq!(
            read_support_index.get_min_read_support(depth),
            min_read_support,
            "depth {}, max FPR {}",
            depth, max_fpr
        );
    }
}


/// At a sequencing error of 0 the null leaves no read to errors, so every depth takes the floor
/// of 2.
#[test]
fn rna_junction_read_support_index_returns_floor_at_zero_sequencing_error() {
    let depths: HashSet<ReadDepth> = HashSet::from([1, 2, 300, 10_000]);

    let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &depths,
        0.0,
        1e-6,
        2
    );

    for depth in depths {
        assert_eq!(read_support_index.get_min_read_support(depth), 2, "depth {}", depth);
    }
}


/// The index is the same on any number of threads, and the same built from a worker of another
/// thread pool, as the splice junction clusterer builds it once per summary group. There it
/// calculates its depths on the worker and builds no pool: a pool built on every call from the
/// workers of another pool exhausts the threads of the process.
///
/// At a sequencing error of 0.01 and a max FPR of 1e-6:
///
///   depth        0    1    2   30   100   1,000   10,000
///   minimum      2    2    3    9    19     138    1,319
#[test]
fn rna_junction_read_support_index_returns_matches_on_any_thread_pool() {
    let depths: HashSet<ReadDepth> = HashSet::from([0, 1, 2, 30, 100, 1_000, 10_000]);

    let on_one_thread: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &depths,
        0.01,
        1e-6,
        1
    );
    let on_four_threads: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
        &depths,
        0.01,
        1e-6,
        4
    );
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let on_worker: Vec<RNAJunctionReadSupportIndex> = thread_pool.install(|| {
        (0..16)
            .into_par_iter()
            .map(|_| RNAJunctionReadSupportIndex::new(
                &depths,
                0.01,
                1e-6,
                4
            ))
            .collect()
    });

    for (depth, min_read_support) in [(0u32, 2u32), (1, 2), (2, 3), (30, 9), (100, 19), (1_000, 138), (10_000, 1_319)] {
        assert_eq!(on_one_thread.get_min_read_support(depth), min_read_support, "depth {}", depth);
    }
    assert_eq!(on_four_threads, on_one_thread);
    for index in on_worker.iter() {
        assert_eq!(index, &on_one_thread);
    }
}
