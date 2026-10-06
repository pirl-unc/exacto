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
#[test]
#[should_panic(expected = "max_k must be > 0")]
fn test_perform_minimum_error_correction_max_k_0() {
    let matrix: ReadMatrix = vec![vec![Some(1.0), Some(0.0)], vec![Some(0.0), Some(1.0)]];
    perform_minimum_error_correction(&matrix, 0, 10, 100, 0);
}

#[test]
#[should_panic(expected = "num_restarts must be > 0")]
fn test_perform_minimum_error_correction_num_restarts_0() {
    let matrix: ReadMatrix = vec![vec![Some(1.0), Some(0.0)], vec![Some(0.0), Some(1.0)]];
    perform_minimum_error_correction(&matrix, 2, 0, 100, 0);
}
