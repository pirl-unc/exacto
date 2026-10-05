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


use crate::stitching::end_resolution::EndResolution;


#[derive(Clone,Debug)]
pub struct StitchedTranscript {
    pub read_id: usize,
    pub original_length: u32,
    pub stitched_length: u32,
    pub stitched_sequence: Box<str>,

    /// `stitched_sequence[..five_prime_stitch_end_index]` is the 5' stitch sequence.
    /// `None` when the 5' end is unresolved.
    pub five_prime_stitch_end_index: Option<usize>,

    /// `stitched_sequence[three_prime_stitch_start_index..]` is the 3' stitch sequence.
    /// `None` when the 3' end is unresolved.
    pub three_prime_stitch_start_index: Option<usize>,

    pub five_prime_end_resolution: Option<EndResolution>,
    pub three_prime_end_resolution: Option<EndResolution>,

    /// Why the read has no transcript model: "mapping_quality" (under the minimum) or
    /// "modeling_error". Such a read passes through as sequenced. `None` when modeled.
    pub unmodeled_reason: Option<&'static str>
}

impl StitchedTranscript {
    pub fn is_five_prime_stitched(&self) -> bool {
        self.five_prime_end_resolution.is_some()
    }

    pub fn is_three_prime_stitched(&self) -> bool {
        self.three_prime_end_resolution.is_some()
    }
}
