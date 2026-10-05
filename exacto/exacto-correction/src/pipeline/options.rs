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


#[derive(Debug, Clone)]
pub struct CorrectRNAReadsOptions {
    pub max_correctable_event_len: usize,
    pub corrected_base_quality: u8,
    /// Cut a read end that runs past an internal annotated exon back to the exon. Off by
    /// default: such an end is also what the clusterer counts as a retained intron.
    pub trim_transcript_ends: bool
}

/// CorrectRNAReadsOptions
impl CorrectRNAReadsOptions {
    pub const DEFAULT: Self = Self {
        max_correctable_event_len: 50,
        corrected_base_quality: 60,
        trim_transcript_ends: false
    };
}
impl Default for CorrectRNAReadsOptions {
    fn default() -> Self {
        Self::DEFAULT
    }
}