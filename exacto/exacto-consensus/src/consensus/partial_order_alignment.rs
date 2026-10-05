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


use abpoa_rs::{
    AlignmentMode,
    AlignmentParameters,
    AlignmentParametersBuilder,
    ConsensusAlgorithm,
    ConsensusData,
    Graph,
    Verbosity
};


pub fn perform_partial_order_alignment(
    reads: &Vec<Box<str>>,
    alignment_mode: AlignmentMode,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32
) -> Box<str> {
    assert!(gap_open_score >= 0, "The gap open score must be at least 0, got {}.", gap_open_score);
    assert!(gap_extend_score >= 1, "The gap extend score must be at least 1, got {}.", gap_extend_score);

    let normalized: Vec<Vec<u8>> = reads
        .iter()
        .filter(|r| !r.is_empty())
        .map(|r| {
            r.bytes()
                .map(|b| match b {
                    b'U' => b'T',
                    b'u' => b't',
                    other => other,
                })
                .collect()
        })
        .collect();

    let alignment_parameters: AlignmentParameters = AlignmentParametersBuilder::new()
        .alignment_mode(alignment_mode)
        .gap_affine_penalties(
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score
        )
        .verbosity(Verbosity::None)
        .build();

    let mut graph = Graph::new(&alignment_parameters);

    let weights: Vec<Vec<i32>> = normalized
        .iter()
        .map(|seq| vec![1; seq.len()]) // uniform per-base weight
        .collect();

    let names: Vec<Vec<u8>> = (0..normalized.len())
        .map(|i| format!("read{}", i + 1).into_bytes())
        .collect();

    graph
        .align_and_add_multiple(
            &alignment_parameters,
            &normalized,
            &weights,
            &names
        )
        .unwrap();

    // One consensus sequence: abPOA reads the minimum frequency only when asked for more.
    graph.generate_consensus_multiple(
        ConsensusAlgorithm::HeaviestBundle,
        1,
        None
    );

    // `get_consensus` is None only if nothing was produced.
    let consensus: ConsensusData = match graph.get_consensus() {
        Some(c) => c,
        None => return "".into()
    };

    let out: String = match consensus.sequences().first() {
        Some(&coded) => String::from_utf8_lossy(
            &alignment_parameters.reverse_seq(coded)
        ).into_owned(),
        None => return "".into()
    };

    out.into()
}


#[cfg(test)]
#[path = "../tests/consensus/partial_order_alignment.rs"]
mod tests;