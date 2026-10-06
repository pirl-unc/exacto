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


use exacto_core::prelude::ReadSupport;

use crate::filtering::rna::template_switch::{
    FoldbackStem, 
    JunctionHomology
};


#[derive(Clone, Debug)]
pub struct TemplateSwitchEvidence {
    /// Panel A: flank homology of a two-sided junction. `None` when the grammar could not
    /// orient the sides or the reference was unavailable.
    pub junction_homology: Option<JunctionHomology>,

    /// Panel B: the inverted repeat at the fold. `None` for every other shape.
    pub foldback_stem: Option<FoldbackStem>,

    /// Both breakends sit on an exon boundary observed in this partition: the shape of a
    /// rearranged or spliced exon, which exempts both panels.
    pub at_exon_boundaries: bool,

    /// Fold-backs only: some supporting read's second arm reaches sequence its first arm
    /// never covered, which a hairpin cannot do.
    pub has_exit_junction: bool,

    /// Pooled breakend families only: member breakpoint spread past the homology interval.
    pub dispersion_beyond_homology: Option<u32>,

    pub num_reads: ReadSupport
}
