use crate::structs::union_find::UnionFind;


#[test]
fn test_union_find_1() {
    let mut uf: UnionFind = UnionFind::new();
    uf.union(1,2);
    uf.union(2,3);
    let clusters: Vec<Vec<usize>> = uf.get_clusters();
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
    let clusters: Vec<Vec<usize>> = uf.get_clusters();
    assert_eq!(clusters.len(), 2);
}