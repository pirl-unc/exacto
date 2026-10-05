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



use super::*;


#[test]
fn test_score_polya_empty_sequence() {
    let score: PolyATailScore = score_polya("");

    assert_eq!(score.total_length, 0, "total length");
    assert_eq!(score.num_adenosine, 0, "number of A bases");
    assert_eq!(score.adenosine_fraction, 0.0, "A fraction is defined as 0.0 for an empty sequence");
    assert_eq!(score.max_consecutive_adenosine, 0, "longest A run");
    assert_eq!(score.trailing_adenosine, 0, "trailing A run");
}


#[test]
fn test_score_polya_all_adenosine() {
    let score: PolyATailScore = score_polya("AAAAAAAAAAAA");

    assert_eq!(score.total_length, 12, "total length");
    assert_eq!(score.num_adenosine, 12, "number of A bases");
    assert_eq!(score.adenosine_fraction, 1.0, "A fraction");
    assert_eq!(score.max_consecutive_adenosine, 12, "longest A run");
    assert_eq!(score.trailing_adenosine, 12, "trailing A run");
}


#[test]
fn test_score_polya_no_adenosine() {
    let score: PolyATailScore = score_polya("GCGCTTGC");

    assert_eq!(score.total_length, 8, "total length");
    assert_eq!(score.num_adenosine, 0, "number of A bases");
    assert_eq!(score.adenosine_fraction, 0.0, "A fraction");
    assert_eq!(score.max_consecutive_adenosine, 0, "longest A run");
    assert_eq!(score.trailing_adenosine, 0, "trailing A run");
}


#[test]
fn test_score_polya_internal_run() {
    // G A T T A C A A A A A G
    // The longest run (5) is internal, and the final base is not an A.
    let score: PolyATailScore = score_polya("GATTACAAAAAG");

    assert_eq!(score.total_length, 12, "total length");
    assert_eq!(score.num_adenosine, 7, "number of A bases");
    assert_eq!(score.adenosine_fraction, 7.0 / 12.0, "A fraction");
    assert_eq!(score.max_consecutive_adenosine, 5, "longest A run");
    assert_eq!(score.trailing_adenosine, 0, "trailing A run");
}


#[test]
fn test_score_polya_trailing_run_shorter_than_longest_run() {
    // The longest run (6) is internal; the trailing run (3) is reported separately.
    let score: PolyATailScore = score_polya("AAAAAAGCAAA");

    assert_eq!(score.total_length, 11, "total length");
    assert_eq!(score.num_adenosine, 9, "number of A bases");
    assert_eq!(score.max_consecutive_adenosine, 6, "longest A run");
    assert_eq!(score.trailing_adenosine, 3, "trailing A run");
}


#[test]
fn test_score_polya_interrupted_tail() {
    // A tail with two sequencing errors still clears a 0.8 A fraction.
    let score: PolyATailScore = score_polya("AAAAAGAAAAAAAACAAAA");

    assert_eq!(score.total_length, 19, "total length");
    assert_eq!(score.num_adenosine, 17, "number of A bases");
    assert!(score.adenosine_fraction >= 0.8, "A fraction {} should clear 0.8", score.adenosine_fraction);
    assert_eq!(score.max_consecutive_adenosine, 8, "longest A run");
    assert_eq!(score.trailing_adenosine, 4, "trailing A run");
}


#[test]
fn test_score_polya_non_adenosine_fragment_before_tail() {
    // A non-A fragment followed by a tail does not clear a 0.8 A fraction.
    let score: PolyATailScore = score_polya("GTCCGTAAAAAAAAAA");

    assert_eq!(score.num_adenosine, 10, "number of A bases");
    assert_eq!(score.adenosine_fraction, 10.0 / 16.0, "A fraction");
    assert!(score.adenosine_fraction < 0.8, "A fraction {} should not clear 0.8", score.adenosine_fraction);
    assert_eq!(score.trailing_adenosine, 10, "trailing A run");
}


#[test]
fn test_score_polya_counts_lowercase() {
    // Soft-masked reference sequence is lowercase.
    let score: PolyATailScore = score_polya("gcaaaAAAgc");

    assert_eq!(score.num_adenosine, 6, "number of A bases");
    assert_eq!(score.max_consecutive_adenosine, 6, "longest A run");
    assert_eq!(score.trailing_adenosine, 0, "trailing A run");
}


#[test]
#[should_panic(expected = "Expected an ASCII nucleotide sequence.")]
fn test_score_polya_non_ascii_sequence_panics() {
    score_polya("AAÅA");
}
