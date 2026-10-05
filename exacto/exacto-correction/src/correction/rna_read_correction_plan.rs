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


use exacto_caller::prelude::{AlignmentModel, AlignmentModelBase};
use std::ops::Range;

use crate::common::enums::BaseAction;


#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BaseCall {
    pub nucleotide: u8,
    pub quality: u8
}


#[derive(Clone, Debug)]
pub(crate) struct BaseEmission {
    /// Emit these bases immediately before this original read position.
    pub before: Vec<BaseCall>,
    pub action: BaseAction
}


#[derive(Clone, Debug)]
pub(crate) struct RNAReadCorrectionPlan {
    /// Exactly one row per original read base.
    pub bases: Vec<BaseEmission>,

    /// Bases emitted after the final original read position.
    pub trailing: Vec<BaseCall>
}

impl RNAReadCorrectionPlan {
    pub fn new(model: &AlignmentModel) -> Self {
        Self {
            bases: vec![
                BaseEmission {
                    before: Vec::new(),
                    action: BaseAction::Keep
                };
                model.get_bases().len()
            ],
            trailing: Vec::new()
        }
    }

    pub fn skip(&mut self, range: Range<usize>) {
        for row in &mut self.bases[range] {
            row.action = BaseAction::Skip;
        }
    }

    pub fn insert_before(
        &mut self,
        boundary: usize,
        calls: impl IntoIterator<Item = BaseCall>,
    ) {
        if boundary == self.bases.len() {
            self.trailing.extend(calls);
        } else {
            self.bases[boundary].before.extend(calls);
        }
    }

    pub fn emit(&self, model: &AlignmentModel) -> (Vec<u8>, Vec<u8>) {
        assert_eq!(self.bases.len(), model.get_bases().len());

        let mut sequence = Vec::new();
        let mut qualities = Vec::new();

        let mut append = |call: BaseCall| {
            sequence.push(call.nucleotide);
            qualities.push(call.quality);
        };

        for (row, base) in self.bases.iter().zip(model.get_bases()) {
            // Boundary additions survive even when this original base is skipped.
            for &call in &row.before {
                append(call);
            }

            match row.action {
                BaseAction::Keep => append(observed_call(base)),
                BaseAction::Replace(call) => append(call),
                BaseAction::Skip => {}
            }
        }

        for &call in &self.trailing {
            append(call);
        }

        (sequence, qualities)
    }
}

pub(crate) fn observed_call(base: &AlignmentModelBase) -> BaseCall {
    BaseCall {
        nucleotide: base.get_nucleotide().as_str().as_bytes()[0],
        quality: base.get_base_quality(),
    }
}