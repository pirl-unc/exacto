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


/// TranscriptModelTopology
#[repr(u8)]
#[derive(Clone,Copy,Debug,Eq,Hash,PartialEq)]
pub enum TranscriptModelTopology {
    /// No fusion or backsplicing.
    Linear,

    /// At least one event with `context = fusion` - a supplementary-alignment
    /// breakpoint or a read-through.
    Fusion,

    /// At least one event row with `context = backsplicing`. Wins over `Fusion` when
    /// both are present: a circular molecule has no free termini to stitch.
    BackSplicing,

    /// The transcript never aligned, so no structure rows exist.
    None
}
impl TranscriptModelTopology {
    pub fn as_str(&self) -> &str {
        match self {
            TranscriptModelTopology::Linear => "linear",
            TranscriptModelTopology::Fusion => "fusion",
            TranscriptModelTopology::BackSplicing => "backsplicing",
            TranscriptModelTopology::None => "none"
        }
    }
}


/// TranscriptTerminus
#[repr(u8)]
#[derive(Clone,Copy,Debug,Eq,Hash,PartialEq)]
pub enum TranscriptTerminusContext {
    Exonic,
    Intronic,
    Intergenic
}
