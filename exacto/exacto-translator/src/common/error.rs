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


use thiserror::Error;


/// Why a translation stopped. Each variant names the transcript, or the file, it stopped on.
#[derive(Debug, Error)]
pub enum TranslatorError {
    #[error("{file} could not be read or written: {reason}")]
    File { file: Box<str>, reason: Box<str> },

    #[error("record {record} of {file} could not be read: {reason}")]
    Record { file: Box<str>, record: usize, reason: Box<str> },

    #[error("{num_unmatched} of {num_transcripts} assembled transcript(s) in the alignments have no support row, \
             starting with: {names}. The support/consensus file and the alignments file describe \
             different transcripts — check that both came from the same run, and that a \
             consensus TSV (transcripts named by cluster id) is not being paired with \
             alignments of externally assembled transcripts (named cid_*/gid_*).")]
    UnmatchedTranscripts { num_unmatched: usize, num_transcripts: usize, names: Box<str> },

    #[error("transcript {transcript}: {field} {value:?} is not a known value")]
    UnknownValue { transcript: Box<str>, field: &'static str, value: Box<str> },

    #[error("transcript {transcript}: its sequence holds a character outside ASCII")]
    NonAsciiSequence { transcript: Box<str> },

    #[error("transcript {transcript}: alignment row {index} (read {read_start}-{read_end}) does not spell \
             the transcript sequence there. The support/consensus file and the alignments file \
             come from different runs.")]
    SequenceMismatch { transcript: Box<str>, index: u32, read_start: u32, read_end: u32 },

    #[error("transcript {transcript}: integrated DNA variant {dna_variant_id} is not in the DNA \
             variants file. The integration and the DNA variants files come from different runs.")]
    UnknownDnaVariant { transcript: Box<str>, dna_variant_id: u32 },

    #[error("transcript {transcript}: sequence length {sequence_length} but its stitch annotation \
             says stitched_length {stitched_length}. Pass the stitched table written by the same \
             stitch-reference-rnas run that produced the stitched sequences.")]
    StitchedLengthMismatch { transcript: Box<str>, sequence_length: usize, stitched_length: u32 },

    #[error("transcript {transcript}: stitch points must satisfy five_prime ({five_prime}) <= \
             three_prime ({three_prime}) <= stitched_length ({stitched_length})")]
    InvalidStitchPoints { transcript: Box<str>, five_prime: u32, three_prime: u32, stitched_length: u32 },

    #[error("the thread pool could not be built: {reason}")]
    ThreadPool { reason: Box<str> }
}
