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


use exacto_core::prelude::*;
use serde::{Deserialize, Serialize};

use crate::prelude::*;
use crate::transcript::transcript_model_annotation::SkippedReferenceRun;


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TranscriptModelEventAnnotation {
    pub(super) context: Option<AlignmentModelEventContext>,

    /// Sorted by (chromosome ID, start); runs never overlap.
    skipped: Vec<SkippedReferenceRun>
}

impl TranscriptModelEventAnnotation {
    pub(super) fn new() -> Self {
        Self {
            context: None,
            skipped: Vec::new()
        }
    }

    pub fn get_context(&self) -> &Option<AlignmentModelEventContext> {
        &self.context
    }
    
    pub fn get_skipped_runs(&self) -> &[SkippedReferenceRun] {
        &self.skipped
    }
    
    pub fn get_skipped_reference_bases(&self) -> Vec<Vec<ReferenceBase>> {
        let mut clusters: Vec<Vec<ReferenceBase>> = Vec::new();
        let mut previous: Option<(ReferenceChromosomeID, ReferencePosition)> = None;
        for run in self.skipped.iter() {
            let adjacent: bool = previous == Some((run.reference_chromosome_id, run.reference_start.wrapping_sub(1)));
            if adjacent {
                clusters.last_mut().unwrap().extend(run.get_bases());
            } else {
                clusters.push(run.get_bases());
            }
            previous = Some((run.reference_chromosome_id, run.reference_end));
        }
        clusters
    }

    pub(super) fn add_skipped_reference_base(&mut self, reference_base: ReferenceBase) {
        let run: SkippedReferenceRun = SkippedReferenceRun::from_reference_base(&reference_base);
        let i: usize = self.skipped.partition_point(|existing| existing.key() < run.key());
        if i > 0 && self.skipped[i - 1].contains(run.reference_chromosome_id, run.reference_start, &run.reference_strand) {
            return;   // already recorded
        }
        self.add_skipped_reference_run(run);
    }

    pub(super) fn add_skipped_reference_run(&mut self, run: SkippedReferenceRun) {
        let i: usize = self.skipped.partition_point(|existing| existing.key() < run.key());
        self.skipped.insert(i, run);
        // Coalesce with the run after, then the run before.
        if i + 1 < self.skipped.len() && self.skipped[i].continues_into(&self.skipped[i + 1]) {
            let next: SkippedReferenceRun = self.skipped.remove(i + 1);
            self.skipped[i].absorb(next);
        }
        if i > 0 && self.skipped[i - 1].continues_into(&self.skipped[i]) {
            let current: SkippedReferenceRun = self.skipped.remove(i);
            self.skipped[i - 1].absorb(current);
        }
    }
}