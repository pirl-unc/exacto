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
fn test_find_polyadenylation_signals_reports_start_and_offset() {
    // sequence    G G G [A A T A A A] C  C  C  C  C  C
    // index       0 1 2  3 4 5 6 7 8  9 10 11 12 13 14
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "GGGAATAAACCCCCC",
        &["AATAAA", "ATTAAA"],
        &(10..=40)
    );

    assert_eq!(hits.len(), 1, "number of hits");
    assert_eq!(&*hits[0].hexamer, "AATAAA", "hexamer");
    assert_eq!(hits[0].read_start, 3, "start of the hexamer in the sequence");

    // Six hexamer bases plus the six bases after them.
    assert_eq!(hits[0].read_end_offset, 12, "distance from the hexamer's start to the sequence's end");
}


#[test]
fn test_find_polyadenylation_signals_no_hexamer_in_sequence() {
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
        &["AATAAA", "ATTAAA"],
        &(10..=40)
    );

    assert!(hits.is_empty(), "a GC-only sequence holds no hexamer");
}


#[test]
fn test_find_polyadenylation_signals_lower_bound_of_range() {
    // The hexamer starts 8 bases before the sequence's end: too close to a read's
    // 3' end to have been used by the cleavage machinery.
    let too_close: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "GCGCAATAAAGC",
        &["AATAAA"],
        &(10..=40)
    );
    assert!(too_close.is_empty(), "an offset of 8 is below the range");

    // The bounds are inclusive.
    let on_bound: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "GCGCAATAAAGC",
        &["AATAAA"],
        &(8..=40)
    );
    assert_eq!(on_bound.len(), 1, "an offset of 8 is on the lower bound");
    assert_eq!(on_bound[0].read_start, 4, "start of the hexamer in the sequence");
    assert_eq!(on_bound[0].read_end_offset, 8, "offset");
}


#[test]
fn test_find_polyadenylation_signals_upper_bound_of_range() {
    // AATAAA followed by 34 bases: offset 40, on the upper bound.
    let on_bound: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "AATAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
        &["AATAAA"],
        &(10..=40)
    );
    assert_eq!(on_bound.len(), 1, "an offset of 40 is on the upper bound");
    assert_eq!(on_bound[0].read_start, 0, "start of the hexamer in the sequence");
    assert_eq!(on_bound[0].read_end_offset, 40, "offset");

    // One more base after it: offset 41, past the upper bound.
    let too_far: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "AATAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        &["AATAAA"],
        &(10..=40)
    );
    assert!(too_far.is_empty(), "an offset of 41 is above the range");
}


#[test]
fn test_find_polyadenylation_signals_nearest_to_the_end_first() {
    // sequence    [A T T A A A] C x10 [A A T A A A] C x10
    // index        0 --------5  6--15  16------21   22--31
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "ATTAAACCCCCCCCCCAATAAACCCCCCCCCC",
        &["AATAAA", "ATTAAA"],
        &(10..=40)
    );

    assert_eq!(hits.len(), 2, "number of hits");
    assert_eq!(&*hits[0].hexamer, "AATAAA", "nearest hexamer");
    assert_eq!(hits[0].read_start, 16, "nearest hexamer start");
    assert_eq!(hits[0].read_end_offset, 16, "nearest hexamer offset");
    assert_eq!(&*hits[1].hexamer, "ATTAAA", "farthest hexamer");
    assert_eq!(hits[1].read_start, 0, "farthest hexamer start");
    assert_eq!(hits[1].read_end_offset, 32, "farthest hexamer offset");
}


#[test]
fn test_find_polyadenylation_signals_only_searches_supplied_hexamers() {
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "ATTAAACCCCCCCCCCAATAAACCCCCCCCCC",
        &["ATTAAA"],
        &(10..=40)
    );

    assert_eq!(hits.len(), 1, "number of hits");
    assert_eq!(&*hits[0].hexamer, "ATTAAA", "hexamer");
    assert_eq!(hits[0].read_end_offset, 32, "offset");

    let none: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "ATTAAACCCCCCCCCCAATAAACCCCCCCCCC",
        &[],
        &(10..=40)
    );
    assert!(none.is_empty(), "no hexamers to search for");
}


#[test]
fn test_find_polyadenylation_signals_normalizes_case_and_uracil() {
    // An RNA spelling in lowercase, searched with an RNA hexamer.
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "gggaauaaacccccc",
        &["aauaaa"],
        &(10..=40)
    );

    assert_eq!(hits.len(), 1, "number of hits");
    assert_eq!(&*hits[0].hexamer, "AATAAA", "the hexamer is reported in uppercase DNA spelling");
    assert_eq!(hits[0].read_start, 3, "start of the hexamer in the sequence");
    assert_eq!(hits[0].read_end_offset, 12, "offset");
}


#[test]
fn test_find_polyadenylation_signals_sequence_shorter_than_a_hexamer() {
    let hits: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "AATAA",
        &["AATAAA"],
        &(1..=40)
    );
    assert!(hits.is_empty(), "five bases cannot hold a hexamer");

    let empty: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        "",
        &["AATAAA"],
        &(1..=40)
    );
    assert!(empty.is_empty(), "an empty sequence holds no hexamer");
}


#[test]
#[should_panic(expected = "PAS motifs must be hexamers.")]
fn test_find_polyadenylation_signals_motif_of_wrong_length_panics() {
    find_polyadenylation_signals("GGGAATAAACCCCCC", &["AATAA"], &(10..=40));
}


#[test]
#[should_panic(expected = "PAS motifs must contain only A/C/G/T/U.")]
fn test_find_polyadenylation_signals_motif_with_ambiguous_base_panics() {
    find_polyadenylation_signals("GGGAATAAACCCCCC", &["AATNAA"], &(10..=40));
}


#[test]
#[should_panic(expected = "Expected ASCII nucleotides.")]
fn test_find_polyadenylation_signals_non_ascii_sequence_panics() {
    find_polyadenylation_signals("GGGAATÅAACCCCCC", &["AATAAA"], &(10..=40));
}
