use super::*;


/// With no shared reads the EM is trivial, so each cluster's CPM is exactly proportional
/// to its read support: 3/8, 3/8, 2/8 of one million.
#[test]
fn quantify_unique_reads_gives_proportional_cpm() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(10, ["a", "b", "c"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(20, ["d", "e", "f"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(30, ["g", "h"].iter().map(|name| Box::<str>::from(*name)).collect());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    let cpm: HashMap<usize, f64> = set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_cpm()))
        .collect();
    assert_eq!(cpm.len(), 3);
    assert!((cpm[&10] - 375_000.0).abs() < 1e-3, "cluster 10 CPM was {}", cpm[&10]);
    assert!((cpm[&20] - 375_000.0).abs() < 1e-3, "cluster 20 CPM was {}", cpm[&20]);
    assert!((cpm[&30] - 250_000.0).abs() < 1e-3, "cluster 30 CPM was {}", cpm[&30]);
}


/// Every input cluster comes back exactly once, carrying its own cluster id and its own
/// read names (not another cluster's), so downstream record builders can join on the id.
#[test]
fn quantify_returns_one_quantification_per_cluster_with_its_reads() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(3, ["a", "b", "c"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(101, ["d"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(7, ["e", "f"].iter().map(|name| Box::<str>::from(*name)).collect());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    assert_eq!(set.quantifications.len(), 3);

    let mut cluster_ids: Vec<usize> = set.quantifications
        .iter()
        .map(|quantification| quantification.get_cluster_id())
        .collect();
    cluster_ids.sort_unstable();
    assert_eq!(cluster_ids, vec![3, 7, 101]);

    for quantification in set.quantifications.iter() {
        let expected: &HashSet<Box<str>> = &clusters[&quantification.get_cluster_id()];
        assert_eq!(quantification.get_read_names(), expected);
        assert_eq!(quantification.get_read_count(), expected.len());
    }
}


/// "amb" is compatible with clusters 10 and 20, which are otherwise symmetric (3 unique
/// reads each), so it splits 0.5/0.5: alpha = 3.5, 3.5, 2 over N = 9 distinct reads. The
/// ambiguous read is still reported in the read names of both clusters it belongs to.
#[test]
fn quantify_splits_ambiguous_read_between_symmetric_clusters() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(10, ["a", "b", "c", "amb"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(20, ["d", "e", "f", "amb"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(30, ["g", "h"].iter().map(|name| Box::<str>::from(*name)).collect());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    let cpm: HashMap<usize, f64> = set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_cpm()))
        .collect();
    assert!((cpm[&10] - cpm[&20]).abs() < 1e-3, "symmetric clusters must get equal CPM");
    assert!((cpm[&10] - 3.5 / 9.0 * 1e6).abs() < 1.0, "cluster 10 CPM was {}", cpm[&10]);
    assert!((cpm[&30] - 2.0 / 9.0 * 1e6).abs() < 1.0, "cluster 30 CPM was {}", cpm[&30]);

    // The shared read is counted once in N, so CPM is conserved.
    let total: f64 = cpm.values().sum();
    assert!((total - 1_000_000.0).abs() < 1e-3, "CPM must sum to 1e6, got {total}");

    // The read-name sets are the input sets, so the ambiguous read shows up in both.
    for quantification in set.quantifications.iter() {
        if quantification.get_cluster_id() == 30 {
            assert_eq!(quantification.get_read_count(), 2);
            assert!(!quantification.get_read_names().contains("amb"));
        } else {
            assert_eq!(quantification.get_read_count(), 4);
            assert!(quantification.get_read_names().contains("amb"));
        }
    }
}


/// Cluster 2's reads are all shared with cluster 1, which also has unique support. Each EM
/// iteration halves cluster 2's share (theta_2 <- 0.5 * theta_2), so the cluster with unique
/// evidence absorbs the ambiguous reads and the nested cluster's CPM goes to ~0, while its
/// read count still reports the two reads that were assigned to it.
#[test]
fn quantify_cluster_without_unique_reads_loses_cpm_to_cluster_with_unique_reads() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(1, ["a", "b", "c", "d"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(2, ["a", "b"].iter().map(|name| Box::<str>::from(*name)).collect());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    let cpm: HashMap<usize, f64> = set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_cpm()))
        .collect();
    assert!(cpm[&1] > 999_999.0, "cluster 1 CPM was {}", cpm[&1]);
    assert!(cpm[&2] >= 0.0 && cpm[&2] < 1.0, "cluster 2 CPM was {}", cpm[&2]);

    let read_counts: HashMap<usize, usize> = set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_read_count()))
        .collect();
    assert_eq!(read_counts[&1], 4);
    assert_eq!(read_counts[&2], 2);
}


/// A single cluster owns every read, so it takes the entire 1e6 of CPM.
#[test]
fn quantify_single_cluster_takes_all_cpm() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(7, ["a", "b", "c", "d"].iter().map(|name| Box::<str>::from(*name)).collect());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    assert_eq!(set.quantifications.len(), 1);
    assert_eq!(set.quantifications[0].get_cluster_id(), 7);
    assert_eq!(set.quantifications[0].get_read_count(), 4);
    assert!((set.quantifications[0].get_cpm() - 1_000_000.0).abs() < 1e-6);
}


/// A cluster with no reads is still reported (the caller asked about it), with zero reads
/// and zero CPM: the pseudo-count only seeds theta, the first M-step hands it no counts.
#[test]
fn quantify_cluster_with_no_reads_gets_zero_cpm() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(1, ["a", "b"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(5, HashSet::new());

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    assert_eq!(set.quantifications.len(), 2);
    for quantification in set.quantifications.iter() {
        if quantification.get_cluster_id() == 5 {
            assert_eq!(quantification.get_read_count(), 0);
            assert_eq!(quantification.get_cpm(), 0.0);
        } else {
            assert_eq!(quantification.get_cluster_id(), 1);
            assert_eq!(quantification.get_read_count(), 2);
            assert!((quantification.get_cpm() - 1_000_000.0).abs() < 1e-6);
        }
    }
}


/// No clusters in, no quantifications out (and no panic on the empty EM).
#[test]
fn quantify_empty_clusters_returns_empty_set() {
    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();

    let set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );

    assert!(set.quantifications.is_empty());
}


/// The options reach the EM. With `max_iter = 0` no EM step runs, so the CPM is the warm
/// start theta_c = (unique_c + pseudo_count) / (sum(unique) + k * pseudo_count): with a
/// pseudo-count of 1 that is (3 + 1) / 6 and (1 + 1) / 6, not the 3/4 and 1/4 the data
/// alone would give.
#[test]
fn quantify_threads_options_through_to_em() {
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(1, ["a", "b", "c"].iter().map(|name| Box::<str>::from(*name)).collect());
    clusters.insert(2, ["d"].iter().map(|name| Box::<str>::from(*name)).collect());

    let options: QuantifyRNAAbundancesOptions = QuantifyRNAAbundancesOptions {
        pseudo_count: 1.0,
        max_iter: 0,
        tol: 1e-9
    };
    let set: ClusterQuantificationSet = quantify_rna_abundances(&clusters, &options);

    let cpm: HashMap<usize, f64> = set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_cpm()))
        .collect();
    assert!((cpm[&1] - 4.0 / 6.0 * 1e6).abs() < 1e-3, "cluster 1 CPM was {}", cpm[&1]);
    assert!((cpm[&2] - 2.0 / 6.0 * 1e6).abs() < 1e-3, "cluster 2 CPM was {}", cpm[&2]);

    // Running the EM to convergence recovers the data-only proportions.
    let converged_set: ClusterQuantificationSet = quantify_rna_abundances(
        &clusters,
        &QuantifyRNAAbundancesOptions::default()
    );
    let converged_cpm: HashMap<usize, f64> = converged_set.quantifications
        .iter()
        .map(|quantification| (quantification.get_cluster_id(), quantification.get_cpm()))
        .collect();
    assert!((converged_cpm[&1] - 750_000.0).abs() < 1e-3, "cluster 1 CPM was {}", converged_cpm[&1]);
    assert!((converged_cpm[&2] - 250_000.0).abs() < 1e-3, "cluster 2 CPM was {}", converged_cpm[&2]);
}
