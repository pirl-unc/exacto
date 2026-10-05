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


use std::ops::RangeInclusive;


#[derive(Clone, Debug)]
pub(crate) struct PolyadenylationSignalMatch {
    /// Example:
    ///
    /// ```text
    /// sequence (input)    G G G [AATAAA] C  C  C  C  C  C
    /// index               0 1 2  3----8  9 10 11 12 13 14
    ///
    /// hexamer             AATAAA
    /// read_start          3
    /// read_end_offset     12
    /// ```
    /// Matched hexamer in uppercase DNA spelling, e.g. AATAAA.
    pub hexamer: Box<str>,

    /// Zero-based position of the hexamer's first base within the supplied read sequence.
    pub read_start: usize,

    /// Distance from the hexamer's first base to the boundary immediately after the read's last base.
    /// Includes the six hexamer bases: a hexamer followed by 14 bases has an offset of 20.
    pub read_end_offset: usize
}


/// Finds exact polyadenylation-signal hexamers near a read's 3′ end.
///
/// The supplied sequence must:
/// - Be oriented 5′ → 3′ along the transcript.
/// - Have its poly(A) tail and terminal adapters already removed.
///
/// `start_offset_range` specifies the permitted distance from the
/// hexamer's START to the boundary immediately after the read.
///
/// Returns all hits, nearest to the read's 3′ end first.
/// These are candidate signals; the read end is not assumed to be
/// a confirmed cleavage site.
pub(crate) fn find_polyadenylation_signals(
    read_sequence: &str,
    hexamers: &[&str],
    start_offset_range: &RangeInclusive<usize>
) -> Vec<PolyadenylationSignalMatch> {
    let normalize = |sequence: &str| -> Vec<u8> {
        assert!(sequence.is_ascii(), "Expected ASCII nucleotides.");

        sequence
            .bytes()
            .map(|base| match base.to_ascii_uppercase() {
                b'U' => b'T',
                base => base,
            })
            .collect()
    };

    let sequence = normalize(read_sequence);

    let motifs: Vec<Vec<u8>> = hexamers
        .iter()
        .map(|hexamer| {
            let motif = normalize(hexamer);

            assert_eq!(motif.len(), 6, "PAS motifs must be hexamers.");
            assert!(
                motif.iter().all(|base| {
                    matches!(*base, b'A' | b'C' | b'G' | b'T')
                }),
                "PAS motifs must contain only A/C/G/T/U."
            );

            motif
        })
        .collect();

    let mut hits: Vec<PolyadenylationSignalMatch> = Vec::new();

    for (read_start, window) in sequence.windows(6).enumerate().rev() {
        let read_end_offset = sequence.len() - read_start;

        // Moving upstream only increases the offset.
        if read_end_offset > *start_offset_range.end() {
            break;
        }

        if !start_offset_range.contains(&read_end_offset) {
            continue;
        }

        if motifs.iter().any(|motif| motif.as_slice() == window) {
            hits.push(PolyadenylationSignalMatch {
                hexamer: String::from_utf8(window.to_vec())
                    .unwrap()
                    .into_boxed_str(),
                read_start,
                read_end_offset
            });
        }
    }

    hits
}


#[cfg(test)]
#[path = "../tests/polyadenylation/signal.rs"]
mod tests;