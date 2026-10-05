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


use exacto_core::prelude::reverse_complement;
use std::cmp::Ordering;
use std::collections::HashSet;


fn encode_kmers(sequence: &str, kmer_size: usize) -> Vec<(u64, u64)> {
    let mask: u64 = u64::MAX >> (64 - 2 * kmer_size);
    let mut kmers: Vec<(u64, u64)> = Vec::new();
    let (mut kmer, mut reverse_kmer, mut num_bases): (u64, u64, usize) = (0, 0, 0);
    for base in sequence.bytes() {
        let code: u64 = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => {
                num_bases = 0;
                continue;
            }
        };
        kmer = ((kmer << 2) | code) & mask;
        reverse_kmer = (reverse_kmer >> 2) | ((3 - code) << (2 * (kmer_size - 1)));
        num_bases += 1;
        if num_bases >= kmer_size {
            kmers.push((kmer, reverse_kmer));
        }
    }
    kmers
}


pub fn orient_reads(reads: &mut [Box<str>], kmer_size: usize) -> usize {
    assert!(kmer_size <= 32, "The orientation k-mer size must be at most 32, got {}.", kmer_size);
    if kmer_size == 0 {
        return 0;
    }

    // Step 1. Collect the k-mers of the longest read that tell its orientations apart.
    let longest_read: &str = match reads.iter().max_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b))) {
        Some(read) => read,
        None => return 0
    };
    let (kmers, reverse_kmers): (HashSet<u64>, HashSet<u64>) = encode_kmers(longest_read, kmer_size)
        .into_iter()
        .unzip();
    let kmers: HashSet<u64> = kmers.difference(&reverse_kmers).copied().collect();

    // Step 2. Compare each read with the longest read.
    // Greater: in the orientation of the longest read. Less: in the other. Equal: undetermined.
    let orientations: Vec<Ordering> = reads
        .iter()
        .map(|read| {
            let (mut num_same, mut num_opposite): (usize, usize) = (0, 0);
            for (kmer, reverse_kmer) in encode_kmers(read, kmer_size) {
                num_same += kmers.contains(&kmer) as usize;
                num_opposite += kmers.contains(&reverse_kmer) as usize;
            }
            num_same.cmp(&num_opposite)
        })
        .collect();

    // Step 3. Reverse-complement the reads of the orientation of fewer reads.
    let num_same: usize = orientations.iter().filter(|&&o| o == Ordering::Greater).count();
    let num_opposite: usize = orientations.iter().filter(|&&o| o == Ordering::Less).count();
    let minority: Ordering = if num_opposite > num_same { Ordering::Greater } else { Ordering::Less };
    let mut num_reversed: usize = 0;
    for (read, orientation) in reads.iter_mut().zip(orientations) {
        if orientation == minority {
            *read = reverse_complement(read);
            num_reversed += 1;
        }
    }
    num_reversed
}


#[cfg(test)]
#[path = "../tests/consensus/read_orientation.rs"]
mod tests;