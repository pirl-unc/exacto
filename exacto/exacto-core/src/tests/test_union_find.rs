use std::collections::BTreeSet;

use crate::prelude::*;


#[test]
fn test_union_find_1() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    let clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].contains(&1u32), true);
    assert_eq!(clusters[0].contains(&2u32), true);
    assert_eq!(clusters[0].contains(&3u32), true);
}

#[test]
fn test_union_find_2() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    uf.union(4,5);
    uf.union(5,6);
    let clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 2);
}

#[test]
fn test_union_find_3() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    let clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}

#[test]
fn test_union_find_4() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    uf.union(1,2);
    let clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}

#[test]
fn test_union_find_5() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    uf.union(1,2);
    uf.union(1,3);
    uf.union(1,4);
    uf.union(1,5);
    uf.union(2,6);
    uf.union(2,7);
    uf.union(2,8);
    uf.union(2,9);
    uf.union(2,10);
    uf.union(1,2);
    let clusters: Vec<BTreeSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}
#[test]
fn test_union_find_6() {
    // The same unions in two union-finds give the same clusters in the same order: members
    // ascending, clusters by smallest member.
    let unions: Vec<(u32, u32)> = (0..50u32).map(|i| ((i * 37) % 100, (i * 37 + 10) % 100)).collect();
    let mut uf_1: UnionFind = UnionFind::new();
    let mut uf_2: UnionFind = UnionFind::new();
    for (x, y) in unions.iter() {
        uf_1.union(*x, *y);
    }
    for (x, y) in unions.iter().rev() {
        uf_2.union(*y, *x);
    }
    let clusters_1: Vec<Vec<u32>> = uf_1.get_clusters().into_iter().map(|cluster| cluster.into_iter().collect()).collect();
    let clusters_2: Vec<Vec<u32>> = uf_2.get_clusters().into_iter().map(|cluster| cluster.into_iter().collect()).collect();
    assert_eq!(clusters_1, clusters_2);
    for cluster in clusters_1.iter() {
        assert!(cluster.windows(2).all(|pair| pair[0] < pair[1]));
    }
    let smallest: Vec<u32> = clusters_1.iter().map(|cluster| cluster[0]).collect();
    assert!(smallest.windows(2).all(|pair| pair[0] < pair[1]));
    let ids: BTreeSet<u32> = unions.iter().flat_map(|(x, y)| [*x, *y]).collect();
    assert_eq!(clusters_1.iter().map(|cluster| cluster.len()).sum::<usize>(), ids.len());
}
