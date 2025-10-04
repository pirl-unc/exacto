use std::collections::HashSet;
use crate::prelude::*;


#[test]
fn test_union_find_1() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    let clusters: Vec<HashSet<usize>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].contains(&1usize), true);
    assert_eq!(clusters[0].contains(&2usize), true);
    assert_eq!(clusters[0].contains(&3usize), true);
}

#[test]
fn test_union_find_2() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    uf.union(4,5);
    uf.union(5,6);
    let clusters: Vec<HashSet<usize>> = uf.get_clusters();
    assert_eq!(clusters.len(), 2);
}

#[test]
fn test_union_find_3() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    let clusters: Vec<HashSet<usize>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}

#[test]
fn test_union_find_4() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    uf.union(1,2);
    let clusters: Vec<HashSet<usize>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}