use super::*;


/// A cluster of more reads than `max_reads_per_cluster` gives that many of its read names, in
/// sorted order, and the same seed gives the same ones. A cluster of no more reads, or a
/// `max_reads_per_cluster` of 0, gives all of them.
///
///   Reads   max_reads_per_cluster   Names returned
///   20      5                       5 of the 20, sorted, the same for seed 42 twice
///   20      20                      all 20, sorted
///   20      0                       all 20, sorted
#[test]
fn subsample_read_names_draws_max_reads_per_cluster_names() {
    let read_names: HashSet<Box<str>> = (0..20)
        .map(|i| format!("read-{:02}", i).into_boxed_str())
        .collect();
    let mut all: Vec<&str> = read_names.iter().map(|read_name| &**read_name).collect();
    all.sort_unstable();

    let subsample: Vec<&str> = subsample_read_names(&read_names, 5, 42);

    assert_eq!(subsample.len(), 5);
    assert!(subsample.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(subsample.iter().all(|read_name| read_names.contains(*read_name)));
    assert_eq!(subsample_read_names(&read_names, 5, 42), subsample);
    assert_ne!(
        (0..10).map(|seed| subsample_read_names(&read_names, 5, seed)).collect::<HashSet<Vec<&str>>>().len(),
        1
    );
    assert_eq!(subsample_read_names(&read_names, 20, 42), all);
    assert_eq!(subsample_read_names(&read_names, 0, 42), all);
}
