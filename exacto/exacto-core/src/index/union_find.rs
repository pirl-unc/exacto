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


use std::collections::{HashMap, HashSet};


pub struct UnionFind {
    // key      =   child ID
    // value    =   (parent ID, size of tree rooted at parent)
    parents: HashMap<u32, (u32, u32)>
}

impl UnionFind {
    pub fn new() -> Self {
        UnionFind {
            parents: HashMap::new()
        }
    }

    /// Get clusters in this union-find.
    ///
    /// # Returns:
    ///
    /// * A vector of clusters where each cluster is a vector of IDs.
    pub fn get_clusters(&mut self) -> Vec<HashSet<u32>> {
        let keys: Vec<u32> = self.parents.keys().cloned().collect();
        let mut map: HashMap<u32, HashSet<u32>> = HashMap::new();
        for key in keys {
            let parent = self.find(key);
            map.entry(parent)
                .or_insert_with(HashSet::new)
                .insert(key);
        }
        map.into_values().collect()
    }

    /// Get size (number of children) in a node (helper function).
    fn get_size(&self, x: u32) -> u32 {
        self.parents.get(&x).map(|(_, size)| *size).unwrap_or(0)
    }

    /// Find the root of an ID with path compression.
    ///
    /// # Parameters:
    ///
    /// * `x` is an ID.
    ///
    /// # Returns:
    ///
    /// * The root ID.
    fn find(&mut self, x: u32) -> u32 {
        if !self.parents.contains_key(&x) {
            self.parents.insert(x, (x, 1));
            return x;
        }

        // Path compression
        let parent: u32 = self.parents.get(&x).unwrap().0;
        if parent != x {
            let root: u32 = self.find(parent);
            self.parents.insert(x, (root, self.get_size(x)));
            root
        } else {
            x
        }
    }

    /// Union (merge) two IDs.
    ///
    /// # Parameters:
    ///
    /// * `x` is an ID.
    /// * `y` is an ID.
    pub fn union(&mut self, x: u32, y: u32) {
        let root_x: u32 = self.find(x);
        let root_y: u32 = self.find(y);

        if root_x != root_y {
            let size_x: u32 = self.get_size(root_x);
            let size_y: u32 = self.get_size(root_y);

            if size_x < size_y {
                self.parents.insert(root_x, (root_y, size_x));
                self.parents.insert(root_y, (root_y, size_x + size_y));
            } else {
                self.parents.insert(root_y, (root_x, size_y));
                self.parents.insert(root_x, (root_x, size_x + size_y));
            }
        }
    }
}
