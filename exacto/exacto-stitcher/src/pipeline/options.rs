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


#[derive(Clone,Debug)]
pub struct StitchReferenceTranscriptsOptions {
    /// Minimum mapping quality for an aligned read (TranscriptModel) to be subject to stitching.
    pub min_mapping_quality: u16,

    /// A reference transcript match with fewer splice junction matches than this is treated as no match.
    pub min_num_splice_junction_matches: usize,

    /// Polyadenylation hexamer signal search size.
    pub pas_search_size: usize,

    /// Polyadenylation hexamers signals.
    pub pas_hexamers: Vec<Box<str>>,

    /// Polyadenylation hexamer start offset range.
    pub pas_start_offset_range: RangeInclusive<usize>,

    /// Polyadenylation tail minimum adenosine fraction.
    pub polya_tail_min_adenosine_fraction: f64,

    /// Polyadenylation window size.
    pub polya_window_size: usize,

    /// Polyadenylation minimum consecutive adenosine.
    pub polya_min_consecutive_adenosine: usize
}

impl Default for StitchReferenceTranscriptsOptions {
    fn default() -> Self {
        Self {
            min_mapping_quality: 0,
            min_num_splice_junction_matches: 1,
            pas_search_size: 40,
            pas_hexamers: vec!["AATAAA".into(), "ATTAAA".into()],
            pas_start_offset_range: 10..=40,
            polya_tail_min_adenosine_fraction: 0.8f64,
            polya_window_size: 20,
            polya_min_consecutive_adenosine: 6
        }
    }
}