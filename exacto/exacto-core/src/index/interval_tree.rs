// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use crate::common::utilities::overlaps;


#[derive(Debug, Clone)]
pub struct Interval<T> {
    pub start: isize,
    pub end: isize,
    pub value: T
}

impl<T> Interval<T> {
    pub fn new(start: isize, end: isize, value: T) -> Self {
        Self { start, end, value }
    }

    pub fn overlaps(&self, start: isize, end: isize) -> bool {
        overlaps(self.start, self.end, start, end)
    }
}

#[derive(Debug, Clone)]
pub struct IntervalTreeNode<T: Clone> {
    pub interval: Interval<T>,
    pub max_end: isize,
    pub height: usize,
    pub left: Option<Box<IntervalTreeNode<T>>>,
    pub right: Option<Box<IntervalTreeNode<T>>>,
}

impl<T: Clone> IntervalTreeNode<T> {
    pub fn new(interval: Interval<T>) -> Self {
        let max_end = interval.end;
        IntervalTreeNode {
            interval,
            max_end,
            height: 1,
            left: None,
            right: None,
        }
    }
}

/// An interval tree kept balanced as an AVL tree, so its depth grows as log(n) whatever the
/// order of insertion.
#[derive(Debug)]
pub struct IntervalTree<T: Clone> {
    pub root: Option<Box<IntervalTreeNode<T>>>
}

impl<T: Clone> IntervalTree<T> {
    pub fn new() -> Self {
        IntervalTree { root: None }
    }

    fn count_nodes(node: &Option<Box<IntervalTreeNode<T>>>) -> usize {
        match node {
            None => 0,
            Some(n) => {
                1 + Self::count_nodes(&n.left) + Self::count_nodes(&n.right)
            }
        }
    }

    pub fn get_size(&self) -> usize {
        Self::count_nodes(&self.root)
    }

    pub fn insert(&mut self, interval: Interval<T>) {
        self.root = Self::insert_node(self.root.take(), interval);
    }

    fn insert_node(
        node: Option<Box<IntervalTreeNode<T>>>,
        interval: Interval<T>,
    ) -> Option<Box<IntervalTreeNode<T>>> {
        match node {
            None => Some(Box::new(IntervalTreeNode::new(interval))),
            Some(mut boxed_node) => {
                // An interval that starts where this one does goes right, after it.
                if interval.start < boxed_node.interval.start {
                    boxed_node.left = Self::insert_node(boxed_node.left.take(), interval);
                } else {
                    boxed_node.right = Self::insert_node(boxed_node.right.take(), interval);
                }
                Some(Self::rebalance(boxed_node))
            }
        }
    }

    fn height(node: &Option<Box<IntervalTreeNode<T>>>) -> usize {
        node.as_ref().map_or(0, |n| n.height)
    }

    /// Recomputes a node's height and max_end from its children.
    fn update(node: &mut Box<IntervalTreeNode<T>>) {
        node.height = 1 + Self::height(&node.left).max(Self::height(&node.right));
        let left_max = node.left.as_ref().map_or(isize::MIN, |left| left.max_end);
        let right_max = node.right.as_ref().map_or(isize::MIN, |right| right.max_end);
        node.max_end = node.interval.end.max(left_max).max(right_max);
    }

    fn rotate_left(mut node: Box<IntervalTreeNode<T>>) -> Box<IntervalTreeNode<T>> {
        let mut right = node.right.take().unwrap();
        node.right = right.left.take();
        Self::update(&mut node);
        right.left = Some(node);
        Self::update(&mut right);
        right
    }

    fn rotate_right(mut node: Box<IntervalTreeNode<T>>) -> Box<IntervalTreeNode<T>> {
        let mut left = node.left.take().unwrap();
        node.left = left.right.take();
        Self::update(&mut node);
        left.right = Some(node);
        Self::update(&mut left);
        left
    }

    /// Restores the AVL balance of a node whose subtrees differ in height by at most 2.
    /// Rotations keep the in-order sequence, so the order of equal starts is kept.
    fn rebalance(mut node: Box<IntervalTreeNode<T>>) -> Box<IntervalTreeNode<T>> {
        Self::update(&mut node);
        let left_height = Self::height(&node.left);
        let right_height = Self::height(&node.right);
        if left_height > right_height + 1 {
            let left = node.left.take().unwrap();
            node.left = Some(if Self::height(&left.left) < Self::height(&left.right) {
                Self::rotate_left(left)
            } else {
                left
            });
            Self::rotate_right(node)
        } else if right_height > left_height + 1 {
            let right = node.right.take().unwrap();
            node.right = Some(if Self::height(&right.right) < Self::height(&right.left) {
                Self::rotate_right(right)
            } else {
                right
            });
            Self::rotate_left(node)
        } else {
            node
        }
    }

    /// Returns the values of the intervals that overlap [start, end], both ends included,
    /// in order of interval start, and in order of insertion for equal starts.
    pub fn overlaps(&self, start: isize, end: isize) -> Vec<&T> {
        assert!(start <= end, "The following must be true: start <= end.");
        let mut results = Vec::new();
        Self::query_overlap(&self.root, start, end, &mut results);
        results
    }

    fn query_overlap<'a>(
        node: &'a Option<Box<IntervalTreeNode<T>>>,
        start: isize,
        end: isize,
        results: &mut Vec<&'a T>,
    ) {
        if let Some(n) = node {
            if let Some(left) = &n.left {
                if left.max_end >= start {
                    Self::query_overlap(&n.left, start, end, results);
                }
            }

            if n.interval.overlaps(start, end) {
                results.push( &n.interval.value);
            }

            if n.interval.start <= end {
                Self::query_overlap(&n.right, start, end, results);
            }
        }
    }
}
