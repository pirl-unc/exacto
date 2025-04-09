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


use std::collections::HashMap;


#[derive(Debug, Default)]
pub struct TrieNode {
    pub nucleotide: char,
    pub is_end: bool,
    pub children: HashMap<char, TrieNode>
}

impl TrieNode {
    fn new(nucleotide: char) -> Self {
        TrieNode {
            nucleotide,
            is_end: false,
            children: HashMap::new(),
        }
    }
}

#[derive(Debug)]
pub struct Trie {
    root: TrieNode
}

impl Trie {
    pub fn new() -> Self {
        Trie {
            root: TrieNode::new(' ')
        }
    }

    /// Check if a sequence exists in this trie.
    pub fn exists(&self, sequence: &str) -> bool {
        self.find_node(sequence).is_some()
    }

    /// Find a node that contains the given sequence (helper function).
    fn find_node(&self, sequence: &str) -> Option<&TrieNode> {
        let mut node = &self.root;
        for nucleotide in sequence.chars() {
            match node.children.get(&nucleotide) {
                Some(next_node) => {
                    node = next_node
                },
                None => {
                    return None
                }
            }
        }
        Some(node)
    }

    /// Insert a sequence.
    pub fn insert(&mut self, sequence: &str) {
        let mut node = &mut self.root;
        for nucleotide in sequence.chars() {
            node = node.children.entry(nucleotide).or_insert_with(|| TrieNode::new(nucleotide));
        }
        node.is_end = true;
    }

    /// Search for all possible sequences in the trie that begin with a prefix.
    pub fn search(&self, prefix: &str) -> Vec<String> {
        let mut sequences = Vec::new();
        if let Some(node) = self.find_node(prefix) {
            self.search_helper(prefix, node, String::new(), &mut sequences);
        }
        sequences
    }

    /// Search for sequences that have the given prefix (helper function).
    fn search_helper(
        &self,
        query: &str,
        current_node: &TrieNode,
        current_prefix: String,
        sequences: &mut Vec<String>,
    ) {
        if current_node.is_end {
            sequences.push(format!("{}{}", query, current_prefix));
        }
        for (&nucleotide, node) in &current_node.children {
            let mut new_prefix = current_prefix.clone();
            new_prefix.push(nucleotide);
            self.search_helper(query, node, new_prefix, sequences);
        }
    }
}
