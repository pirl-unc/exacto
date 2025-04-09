use crate::structs::interval_tree::{
    Interval,IntervalTree
};


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