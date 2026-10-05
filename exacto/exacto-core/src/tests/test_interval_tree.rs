use crate::prelude::*;


#[test]
fn test_interval_tree_1() {
    let mut itree: IntervalTree<usize> = IntervalTree::new();
    itree.insert(Interval::new(1,5,100));
    itree.insert(Interval::new(6,10,200));
    itree.insert(Interval::new(11,15,300));
    itree.insert(Interval::new(16,20,400));
    itree.insert(Interval::new(21,25,500));

    assert!(itree.overlaps(1,2).len() == 1);
    assert!(*itree.overlaps(1,2)[0] == 100);
    assert!(itree.overlaps(5,12).len() == 3);
    assert!(itree.overlaps(15,25).len() == 3);
}

#[test]
fn test_interval_tree_2() {
    let mut itree: IntervalTree<usize> = IntervalTree::new();
    itree.insert(Interval::new(100,200,100));
    itree.insert(Interval::new(200,1000,200));
    itree.insert(Interval::new(400,1000,300));
    itree.insert(Interval::new(500,1000,400));
    itree.insert(Interval::new(2000,3000,500));

    assert!(itree.overlaps(50,150).len() == 1);
    assert!(*itree.overlaps(50,150)[0] == 100);
    assert!(itree.overlaps(2000,2000).len() == 1);
    assert!(*itree.overlaps(2000,2000)[0] == 500);
}

#[test]
fn test_interval_tree_3() {
    let mut itree: IntervalTree<&str> = IntervalTree::new();
    itree.insert(Interval::new(100,200,"A"));
    itree.insert(Interval::new(200,1000,"B"));
    itree.insert(Interval::new(400,1000,"C"));
    itree.insert(Interval::new(500,1000,"D"));
    itree.insert(Interval::new(2000,3000,"E"));

    assert!(itree.overlaps(50,150).len() == 1);
    assert!(*itree.overlaps(50,150)[0] == "A");
    assert!(itree.overlaps(2000,2000).len() == 1);
    assert!(*itree.overlaps(2000,2000)[0] == "E");
}

#[test]
fn test_interval_tree_4() {
    #[derive(Debug, Clone)]
    struct Test {
        value: String
    }

    let mut itree: IntervalTree<Test> = IntervalTree::new();
    itree.insert(Interval::new(100,200,Test{ value: "A".to_string() }));
    itree.insert(Interval::new(200,1000,Test{ value: "B".to_string() }));
    itree.insert(Interval::new(400,1000,Test{ value: "C".to_string() }));
    itree.insert(Interval::new(500,1000,Test{ value: "D".to_string() }));
    itree.insert(Interval::new(2000,3000,Test{ value: "E".to_string() }));

    assert!(itree.overlaps(50,150).len() == 1);
    assert!(*itree.overlaps(50,150)[0].value == "A".to_string());
    assert!(itree.overlaps(2000,2000).len() == 1);
    assert!(*itree.overlaps(2000,2000)[0].value == "E".to_string());
}
#[test]
fn test_interval_tree_5() {
    // Intervals inserted in order of start (the order exacto-caller builds its breakpoint trees in)
    // and in a scrambled order give the same answers as a scan of every interval: the intervals
    // that overlap the query, in order of start, and in order of insertion for equal starts.
    let mut intervals: Vec<(isize, isize, usize)> = Vec::new();
    for i in 0..2_000 {
        let start: isize = (i / 3) as isize * 7;
        let end: isize = start + ((i * 37) % 50) as isize;
        intervals.push((start, end, i));
    }
    let mut scrambled: Vec<(isize, isize, usize)> = intervals.clone();
    let mut state: usize = 1;
    for i in (1..scrambled.len()).rev() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        scrambled.swap(i, (state >> 33) % (i + 1));
    }
    for order in [intervals, scrambled] {
        let mut itree: IntervalTree<usize> = IntervalTree::new();
        for (start, end, value) in order.iter() {
            itree.insert(Interval::new(*start, *end, *value));
        }
        assert_eq!(itree.get_size(), 2_000);
        for query_start in (-10..4_700).step_by(13) {
            for query_length in [0, 5, 60] {
                let query_end: isize = query_start + query_length;
                let mut expected: Vec<(isize, usize)> = Vec::new();
                for (start, end, value) in order.iter() {
                    if *start <= query_end && *end >= query_start {
                        expected.push((*start, *value));
                    }
                }
                expected.sort_by_key(|(start, _)| *start);
                let expected: Vec<usize> = expected.into_iter().map(|(_, value)| value).collect();
                let results: Vec<usize> = itree.overlaps(query_start, query_end).into_iter().copied().collect();
                assert_eq!(results, expected);
            }
        }
    }
}

#[test]
fn test_interval_tree_6() {
    // 200,000 intervals in order of start, on a thread with a 256 KiB stack: the tree stays
    // balanced (an AVL tree of n nodes is at most 1.44 log2(n) deep), so neither the build nor
    // the queries recurse deeper than that.
    let handle = std::thread::Builder::new()
        .stack_size(256 * 1024)
        .spawn(|| {
            let mut itree: IntervalTree<usize> = IntervalTree::new();
            for i in 0..200_000 {
                itree.insert(Interval::new(i as isize * 10, i as isize * 10 + 15, i));
            }
            let height: usize = itree.root.as_ref().unwrap().height;
            let results: Vec<usize> = itree.overlaps(1_000_000, 1_000_020).into_iter().copied().collect();
            (itree.get_size(), height, results)
        })
        .unwrap();
    let (size, height, results) = handle.join().unwrap();
    assert_eq!(size, 200_000);
    assert!(height <= 25, "height {height}");
    assert_eq!(results, vec![99_999, 100_000, 100_001, 100_002]);
}
