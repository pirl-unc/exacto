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


#[derive(Clone, Debug)]
pub(crate) struct PolyATailScore {
    /// Total number of bases in the supplied sequence including non-A bases.
    pub total_length: usize,

    /// Number of A bases anywhere in the sequence.
    /// Both uppercase A and lowercase a are counted.
    pub num_adenosine: usize,

    /// Fraction of bases that are A: num_adenosine / length.
    /// Ranges from 0.0 to 1.0; defined as 0.0 for an empty sequence.
    pub adenosine_fraction: f64,

    /// Length of the longest uninterrupted adenosine run anywhere
    /// in the sequence.
    pub max_consecutive_adenosine: usize,

    /// Length of the uninterrupted A run at the sequence's end.
    /// Zero if the final base is not A or the sequence is empty.
    pub trailing_adenosine: usize
}


pub(crate) fn score_polya(sequence: &str) -> PolyATailScore {
    assert!(sequence.is_ascii(), "Expected an ASCII nucleotide sequence.");

    let mut num_adenosine: usize = 0;
    let mut run: usize = 0;
    let mut max_consecutive_adenosine: usize = 0;

    for base in sequence.bytes() {
        if base.to_ascii_uppercase() == b'A' {
            num_adenosine += 1;
            run += 1;
            max_consecutive_adenosine = max_consecutive_adenosine.max(run);
        } else {
            run = 0;
        }
    }

    let length: usize = sequence.len();

    PolyATailScore {
        total_length: length,
        num_adenosine: num_adenosine,
        adenosine_fraction: if length == 0 {
            0.0
        } else {
            num_adenosine as f64 / length as f64
        },
        max_consecutive_adenosine,
        trailing_adenosine: run
    }
}


#[cfg(test)]
#[path = "../tests/polyadenylation/tail.rs"]
mod tests;