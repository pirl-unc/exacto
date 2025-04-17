use std::any::Any;
use crate::prelude::*;


#[test]
fn test_calculate_cosine_similarity_1() {
    let v1: Vec<i8> = vec![1, 0, 1, 0];
    let v2: Vec<i8> = vec![1, 0, 1, 0];
    let cosine_similarity_score: f64 = calculate_cosine_similarity(&v1, &v2);
    assert!(cosine_similarity_score > 0.999f64);
}

#[test]
fn test_calculate_cosine_similarity_2() {
    let v1: Vec<i8> = vec![1, 0, 1, 0];
    let v2: Vec<i8> = vec![0, 1, 0, 1];
    let cosine_similarity_score: f64 = calculate_cosine_similarity(&v1, &v2);
    assert!(cosine_similarity_score < 0.001f64);
}

#[test]
fn test_clone_boxed_any_1() {
    let value1: Box<str> = "helloworld".to_string().into();
    let boxed_value1: Box<dyn Any> = Box::new(value1); // Box the value as Box<dyn Any>
    let value2: Box<dyn Any> = clone_boxed_any(&boxed_value1);
    if let Some(cloned_boxed_str) = value2.downcast_ref::<Box<str>>() {
        assert_eq!(&**cloned_boxed_str, "helloworld");
    } else {
        panic!("Unexpected error.");
    }
}

#[test]
fn test_overlaps_1() {
    assert_eq!(overlaps(100,200,150,250), true);
    assert_eq!(overlaps(100,200,99,100), true);
    assert_eq!(overlaps(100,200,200,201), true);
    assert_eq!(overlaps(1000,2000,3000,4000), false);
    assert_eq!(overlaps(5,6,6,7), true);
    assert_eq!(overlaps(1,2000,500,600), true);
    assert_eq!(overlaps(500,600,1,2000), true);
    assert_eq!(overlaps(500,600,550,2000), true);
    assert_eq!(overlaps(1000,2000,550,1100), true);
}

#[test]
fn test_find_overlap_1() {
    assert!(find_overlap((100,200), (150,250)).is_some());
    assert!(find_overlap((100,200), (150,250)).unwrap().0 == 150);
    assert!(find_overlap((100,200), (150,250)).unwrap().1 == 200);
    assert!(find_overlap((5,6), (6,7)).is_some());
    assert!(find_overlap((5,6), (6,7)).unwrap().0 == 6);
    assert!(find_overlap((5,6), (6,7)).unwrap().1 == 6);
    assert!(find_overlap((1,2000), (500,600)).is_some());
    assert!(find_overlap((1,2000), (500,600)).unwrap().0 == 500);
    assert!(find_overlap((1,2000), (500,600)).unwrap().1 == 600);
    assert!(find_overlap((500,600), (1,2000)).is_some());
    assert!(find_overlap((500,600), (1,2000)).unwrap().0 == 500);
    assert!(find_overlap((500,600), (1,2000)).unwrap().1 == 600);
}

#[test]
fn test_merge_regions_1() {
    let mut regions: Vec<(isize,isize)> = Vec::new();
    regions.push((1,5));
    regions.push((2,6));
    regions.push((8,10));
    regions.push((9,12));
    let merged_regions: Vec<(isize,isize)> = merge_regions(regions);
    assert!(merged_regions.len() == 2);
    assert!(merged_regions[0] == (1,6));
    assert!(merged_regions[1] == (8,12));
}

#[test]
fn test_merge_regions_2() {
    let mut regions: Vec<(isize,isize)> = Vec::new();
    regions.push((1,100));
    regions.push((1,5));
    regions.push((2,200));
    regions.push((199,200));
    regions.push((300,400));
    let merged_regions: Vec<(isize,isize)> = merge_regions(regions);
    assert!(merged_regions.len() == 2);
    assert!(merged_regions[0] == (1,200));
    assert!(merged_regions[1] == (300,400));
}

#[test]
fn test_merge_regions_3() {
    let mut regions: Vec<(isize,isize)> = Vec::new();
    regions.push((1,100));
    regions.push((500,501));
    regions.push((300,500));
    let merged_regions: Vec<(isize,isize)> = merge_regions(regions);
    assert!(merged_regions.len() == 2);
    assert!(merged_regions[0] == (1,100));
    assert!(merged_regions[1] == (300,501));
}

#[test]
fn test_count_common_bases_1() {
    let a = vec![("chr1".into(), 1, 100),("chr1".into(),201,300)];
    let b = vec![("chr1".into(), 1, 50),("chr1".into(),201,250)];
    let unioned_bases: u32 = count_common_bases(&a,&b);
    assert_eq!(unioned_bases, 100);
}

#[test]
fn test_count_common_bases_2() {
    let a = vec![("chr1".into(), 1, 100),("chr1".into(),201,300)];
    let b = vec![("chr1".into(), 201, 250),("chr1".into(),291,300)];
    let num_common_bases: u32 = count_common_bases(&a,&b);
    assert_eq!(num_common_bases, 60);
}

#[test]
fn test_count_union_bases_1() {
    let a = vec![("chr1".into(), 1, 100),("chr1".into(),201,300)];
    let b = vec![("chr1".into(), 250, 300),("chr1".into(),401,500)];
    let num_unioned_bases: u32 = count_union_bases(&a,&b);
    assert_eq!(num_unioned_bases, 300);
}

#[test]
fn test_count_union_bases_2() {
    let a = vec![("chr1".into(), 1, 100),("chr1".into(),201,300)];
    let b = vec![("chr2".into(), 1, 100),("chr1".into(),401,500)];
    let num_unioned_bases: u32 = count_union_bases(&a,&b);
    assert_eq!(num_unioned_bases, 400);
}

#[test]
fn test_count_non_overlapping_bases() {
    let a = vec![
        (Box::from("chr1"), 10, 20),
        (Box::from("chr1"), 30, 40),
    ];

    let b = vec![
        (Box::from("chr1"), 15, 35),
    ];

    let (num_a_only_bases, num_b_only_bases) = count_non_overlapping_bases(&a, &b);
    let num_non_overlapping_bases: u32 = num_a_only_bases + num_b_only_bases;
    assert_eq!(num_non_overlapping_bases, 19);
}