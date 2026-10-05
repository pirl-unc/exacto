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
fn test_calculate_l2_distance_1() {
    let v1: Vec<i8> = vec![1, 0, 1, 0];
    let v2: Vec<i8> = vec![0, 1, 0, 1];
    let l2_distance: f32 = calculate_l2_distance(&v1, &v2);
    assert!(l2_distance == 2f32);
}

#[test]
fn test_calculate_l2_distance_2() {
    let v1: Vec<i8> = vec![1, 0, 1, 0];
    let v2: Vec<i8> = vec![1, 1, 1, 0];
    let l2_distance: f32 = calculate_l2_distance(&v1, &v2);
    assert!(l2_distance == 1f32);
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
fn test_count_non_overlapping_bases_1() {
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
fn test_interval_contains_1() {
    assert_eq!(interval_contains(1, 100, 50, 60), true);
}

#[test]
fn test_interval_contains_2() {
    assert_eq!(interval_contains(1, 100, 50, 100), true);
}

#[test]
fn test_interval_contains_3() {
    assert_eq!(interval_contains(1, 100, 50, 101), false);
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
fn test_count_non_overlapping_bases_2() {
    // Positions covered by neither list, one list, or both, plus abutting and
    // duplicated intervals and a chromosome that appears on one side only.
    let a = vec![
        (Box::from("chr1"), 10u32, 20u32),
        (Box::from("chr1"), 21, 25),   // abuts the interval above
        (Box::from("chr1"), 12, 18),   // wholly inside it
        (Box::from("chr2"), 100, 109),
    ];
    let b = vec![
        (Box::from("chr1"), 20, 22),
        (Box::from("chr3"), 1, 5),
    ];

    // chr1: a covers 10-25 (16), b covers 20-22 (3), shared 3 -> a-only 13, b-only 0
    // chr2: a only, 10 bases
    // chr3: b only, 5 bases
    let (num_a_only_bases, num_b_only_bases) = count_non_overlapping_bases(&a, &b);
    assert_eq!(num_a_only_bases, 23);
    assert_eq!(num_b_only_bases, 5);
}

#[test]
fn test_count_non_overlapping_bases_3() {
    // Randomised agreement with a per-position oracle. The production implementation is
    // interval arithmetic; this is the definition it has to match, checked over a coordinate
    // range small enough to enumerate. Deterministic PRNG so a failure is reproducible.
    use std::collections::HashSet;

    fn oracle(
        a: &Vec<(Box<str>, u32, u32)>,
        b: &Vec<(Box<str>, u32, u32)>
    ) -> (u32, u32) {
        fn positions(regions: &Vec<(Box<str>, u32, u32)>) -> HashSet<(Box<str>, u32)> {
            let mut covered: HashSet<(Box<str>, u32)> = HashSet::new();
            for (chromosome, start, end) in regions.iter() {
                for position in *start..=*end {
                    covered.insert((chromosome.clone(), position));
                }
            }
            covered
        }

        let a_positions: HashSet<(Box<str>, u32)> = positions(a);
        let b_positions: HashSet<(Box<str>, u32)> = positions(b);
        (
            a_positions.difference(&b_positions).count() as u32,
            b_positions.difference(&a_positions).count() as u32
        )
    }

    fn next(seed: &mut u64, bound: u32) -> u32 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 33) as u32) % bound
    }

    fn build(seed: &mut u64) -> Vec<(Box<str>, u32, u32)> {
        let num_regions: u32 = next(seed, 6);
        let mut regions: Vec<(Box<str>, u32, u32)> = Vec::new();
        for _ in 0..num_regions {
            let chromosome: Box<str> = if next(seed, 3) == 0 { "chr2".into() } else { "chr1".into() };
            let start: u32 = next(seed, 60);
            let end: u32 = start + next(seed, 20);
            regions.push((chromosome, start, end));
        }
        regions
    }

    let mut seed: u64 = 0x5DEECE66D;
    for _trial in 0..500 {
        let a: Vec<(Box<str>, u32, u32)> = build(&mut seed);
        let b: Vec<(Box<str>, u32, u32)> = build(&mut seed);

        assert_eq!(
            count_non_overlapping_bases(&a, &b),
            oracle(&a, &b),
            "disagreement on a={:?} b={:?}",
            a, b
        );
    }
}

#[test]
fn test_write_tsv_table_empty() {
    // A table without rows is written as its header line, which any TSV reader reads back as an
    // empty table; write_tsv_file leaves a file of 0 bytes.
    #[derive(serde::Serialize, Default)]
    struct Row {
        variant_id: Box<str>,
        depth: u32
    }
    let temp_dir = tempfile::TempDir::new().unwrap();
    let table_file = temp_dir.path().join("table.tsv");
    write_tsv_table(Vec::<Row>::new(), &table_file).unwrap();
    assert_eq!(std::fs::read_to_string(&table_file).unwrap(), "variant_id\tdepth\n");
    write_tsv_table(vec![Row { variant_id: "v1".into(), depth: 3 }], &table_file).unwrap();
    assert_eq!(std::fs::read_to_string(&table_file).unwrap(), "variant_id\tdepth\nv1\t3\n");
}
