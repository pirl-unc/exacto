use crate::prelude::*;


#[test]
fn test_sweep_overlaps() {
    let segments: Vec<(u32, u32)> = vec![(1, 3), (3, 4), (4,6), (15, 20), (25, 30), (27,30)];
    let pairs = sweep_overlaps(&segments);
    assert!(pairs.len() == 3);
    println!("{:?}", pairs);
    assert!(pairs[0] == (0,1));
    assert!(pairs[1] == (1,2));
    assert!(pairs[2] == (4,5));
}