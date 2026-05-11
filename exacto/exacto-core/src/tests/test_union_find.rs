use std::collections::HashSet;

use crate::prelude::*;


#[test]
fn test_union_find_1() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    let clusters: Vec<HashSet<u32>> = uf.get_clusters();
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
    let clusters: Vec<HashSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 2);
}

#[test]
fn test_union_find_3() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    let clusters: Vec<HashSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}

#[test]
fn test_union_find_4() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,1);
    uf.union(1,2);
    let clusters: Vec<HashSet<u32>> = uf.get_clusters();
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
    let clusters: Vec<HashSet<u32>> = uf.get_clusters();
    assert_eq!(clusters.len(), 1);
}