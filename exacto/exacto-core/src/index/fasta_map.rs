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
use std::sync::Arc;

use crate::common::fasta::*;


pub struct FastaMap {
    contigs: HashMap<Box<str>, Arc<[u8]>>
}

impl FastaMap {
    pub fn new(fasta_file: &str) -> Self {
        let mut contigs: HashMap<Box<str>, Arc<[u8]>> = HashMap::new();
        let sequences: Vec<(Box<str>, u32)> = get_fasta_sequence_ids(fasta_file);
        for (contig, len) in sequences {
            let sequence: Box<str> = get_fasta_sequence(&*contig, 1, len, fasta_file);
            let bytes: Vec<u8> = sequence.as_bytes().to_vec();
            contigs.insert(contig, Arc::from(bytes.into_boxed_slice()));
        }

        Self {
            contigs: contigs
        }
    }

    pub fn get_length(&self, contig: &str) -> usize {
        self.contigs.get(contig).unwrap().len()
    }

    /// Get a 1-based inclusive slice from a contig as &str.
    /// Example: slice_str("chr1", 1, 3) -> bases 1..3 (inclusive).
    pub fn get_sequence(&self, contig: &str, start_1b: usize, end_1b: usize) -> &str {
        let seq = self.contigs.get(contig).unwrap();

        // Convert from 1-based inclusive to 0-based [start, end)
        let start: usize = start_1b.checked_sub(1).unwrap();
        let end: usize = end_1b;

        if start >= seq.len() || end > seq.len() || start > end {
            panic!("Invalid range: [{}, {})", start_1b, end_1b);
        }

        std::str::from_utf8(&seq[start..end]).unwrap()
    }
}
