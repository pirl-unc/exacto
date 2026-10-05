use super::*;


fn read_set(names: &[&str]) -> HashSet<Box<str>> {
    names.iter().map(|name| Box::<str>::from(*name)).collect()
}


fn cpm_by_cluster(cpm: &[(usize, f64)]) -> HashMap<usize, f64> {
    cpm.iter().cloned().collect()
}


/// With no ambiguity the E-step is trivial, so each cluster's CPM is exactly
/// proportional to its (unique) read support: 3/8, 3/8, 2/8 of one million.
#[test]
fn em_all_unique_reads_gives_proportional_cpm() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(10, read_set(&["a", "b", "c"]));
    clusters.insert(20, read_set(&["d", "e", "f"]));
    clusters.insert(30, read_set(&["g", "h"]));

    let (converged, cpm, log_likelihoods) = run_finite_mixture_em(&clusters, 1e-6, 1000, 1e-9);

    assert!(converged, "EM should converge on a well-separated input");

    let cpm = cpm_by_cluster(&cpm);
    assert!((cpm[&10] - 375_000.0).abs() < 1e-3);
    assert!((cpm[&20] - 375_000.0).abs() < 1e-3);
    assert!((cpm[&30] - 250_000.0).abs() < 1e-3);

    // CPM is conserved: the whole set sums to one million.
    let total: f64 = cpm.values().sum();
    assert!((total - 1_000_000.0).abs() < 1e-3, "CPM must sum to 1e6, got {total}");

    // EM never decreases the log-likelihood between iterations.
    assert!(!log_likelihoods.is_empty());
    for window in log_likelihoods.windows(2) {
        assert!(
            window[1] >= window[0] - 1e-9,
            "log-likelihood decreased: {window:?}"
        );
    }
}


/// "amb" is compatible with both cluster 10 and 20, which are otherwise symmetric
/// (3 unique reads each). By symmetry the ambiguous read splits 0.5/0.5, so 10 and
/// 20 receive equal CPM (alpha 3.5 each) and 30 receives less (alpha 2), over N = 9.
#[test]
fn em_splits_ambiguous_read_between_symmetric_clusters() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(10, read_set(&["a", "b", "c", "amb"]));
    clusters.insert(20, read_set(&["d", "e", "f", "amb"]));
    clusters.insert(30, read_set(&["g", "h"]));

    let (converged, cpm, _) = run_finite_mixture_em(&clusters, 1e-6, 1000, 1e-9);
    assert!(converged);

    let cpm = cpm_by_cluster(&cpm);
    // Symmetry: 10 and 20 get the same share.
    assert!((cpm[&10] - cpm[&20]).abs() < 1e-3, "symmetric clusters must be equal");
    // 30 has strictly less support than 10.
    assert!(cpm[&30] < cpm[&10]);
    // Expected shares: alpha = 3.5, 3.5, 2 over N = 9.
    assert!((cpm[&10] - 3.5 / 9.0 * 1e6).abs() < 1.0);
    assert!((cpm[&30] - 2.0 / 9.0 * 1e6).abs() < 1.0);

    let total: f64 = cpm.values().sum();
    assert!((total - 1_000_000.0).abs() < 1e-3);
}


/// A single cluster owns every read, so it takes the entire 1e6 of CPM.
#[test]
fn em_single_cluster_takes_all_cpm() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(7, read_set(&["a", "b", "c", "d"]));

    let (converged, cpm, _) = run_finite_mixture_em(&clusters, 1e-6, 1000, 1e-9);
    assert!(converged);

    let cpm = cpm_by_cluster(&cpm);
    assert!((cpm[&7] - 1_000_000.0).abs() < 1e-6);
}
