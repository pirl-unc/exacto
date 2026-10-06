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
use std::ops::Range;

use crate::prelude::*;


pub struct AlignmentModelEditor<'a> {
    model: &'a mut AlignmentModel,

    /// One flag per base of the model: true when this editor put the base in.
    restored: Vec<bool>,

    /// Edits applied, oldest first.
    history: Vec<AlignmentModelEdit>
}


#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AlignmentModelEdit {
    Substitution { read_position: ReadPosition, nucleotide: Nucleotide, quality: BaseQuality },
    Removal { range: Range<ReadPosition> },
    DeletionRestoration { read_position: ReadPosition, quality: BaseQuality },
    UnplacedInsertion { read_position: ReadPosition, nucleotides: Vec<Nucleotide>, quality: BaseQuality },
    Retention { range: Range<ReadPosition> }
}

impl<'a> AlignmentModelEditor<'a> {
    pub fn new(model: &'a mut AlignmentModel) -> Self {
        let restored: Vec<bool> = vec![false; model.num_bases() as usize];
        Self { model, restored, history: Vec::new() }
    }

    pub fn get_model(&self) -> &AlignmentModel {
        self.model
    }

    pub fn get_history(&self) -> &[AlignmentModelEdit] {
        &self.history
    }

    pub fn is_restored(&self, read_position: ReadPosition) -> bool {
        self.restored[read_position as usize]
    }

    pub fn substitute(&mut self, read_position: ReadPosition, nucleotide: Nucleotide, quality: BaseQuality) -> bool {
        let base: &mut AlignmentModelBase = self.model.get_base_mut(read_position);
        if *base.get_nucleotide() == nucleotide {
            return false;
        }
        base.replace_nucleotide(nucleotide.clone(), quality);
        self.history.push(AlignmentModelEdit::Substitution { read_position, nucleotide, quality });
        true
    }
    
    pub fn remove(&mut self, range: Range<ReadPosition>) {
        self.splice(range.clone(), Vec::new());
        self.history.push(AlignmentModelEdit::Removal { range });
    }
    
    pub fn restore_deletion(&mut self, read_position: ReadPosition, quality: BaseQuality) -> bool {
        let Some(flank) = [read_position.checked_sub(1), Some(read_position)]
            .into_iter()
            .flatten()
            .find(|&flank| {
                self.model
                    .get_bases()
                    .get(flank as usize)
                    .is_some_and(|base| base.get_deletion_read_position() == Some(read_position))
            })
        else {
            return false;
        };
        let base: &AlignmentModelBase = self.model.get_base(flank);
        let Some((chromosome_id, position, strand)) = base.get_placement().get_coordinate() else {
            return false;
        };
        let strand: Strand = strand.clone();
        let mapping_quality: MappingQuality = base.get_placement().get_mapping_quality().unwrap_or(0);
        let deleted: Vec<Nucleotide> = base.get_deleted_reference_bases().to_vec();
        let length: u32 = deleted.len() as u32;
        let restored: Vec<AlignmentModelBase> = deleted
            .into_iter()
            .enumerate()
            .map(|(index, nucleotide)| {
                let index: u32 = index as u32;
                let offset: u32 = if strand == Strand::Forward { index + 1 } else { length - index };
                let mut base = AlignmentModelBase::new(read_position + index, nucleotide.clone(), quality);
                base.set_kind(AlignmentModelBaseKind::Match);
                base.place(chromosome_id, position + offset, strand.clone(), mapping_quality);
                base.set_reference_nucleotide(Some(nucleotide));
                base
            })
            .collect();

        self.model.get_base_mut(flank).set_deleted_reference_bases(Vec::new());
        let forward: bool = strand == Strand::Forward;
        self.model.retain_events(|event| {
            *event.get_kind() != AlignmentModelEventKind::Deletion
                || if forward {
                event.get_prev_read_position() != flank
            } else {
                event.get_next_read_position() != flank
            }
        });
        self.splice(read_position..read_position, restored);
        self.history.push(AlignmentModelEdit::DeletionRestoration { read_position, quality });
        true
    }

    pub fn insert_unplaced(&mut self, read_position: ReadPosition, nucleotides: &[Nucleotide], quality: BaseQuality) {
        let inserted: Vec<AlignmentModelBase> = nucleotides
            .iter()
            .enumerate()
            .map(|(index, nucleotide)| {
                let mut base = AlignmentModelBase::new(read_position + index as u32, nucleotide.clone(), quality);
                base.set_kind(AlignmentModelBaseKind::Softclip);
                base
            })
            .collect();
        self.splice(read_position..read_position, inserted);
        self.history.push(AlignmentModelEdit::UnplacedInsertion { read_position, nucleotides: nucleotides.to_vec(), quality });
    }
    
    pub fn retain(&mut self, range: Range<ReadPosition>) {
        let length: u32 = self.model.num_bases();
        assert!(range.start <= range.end && range.end <= length);
        for base in self.model.get_bases_mut() {
            if base.get_deletion_read_position().is_some_and(|position| {
                (range.start > 0 && position <= range.start) || (range.end < length && position >= range.end)
            }) {
                base.set_deleted_reference_bases(Vec::new());
            }
        }
        self.splice(range.end..length, Vec::new());
        self.splice(0..range.start, Vec::new());
        self.history.push(AlignmentModelEdit::Retention { range });
    }

    fn splice(&mut self, range: Range<ReadPosition>, replacement: Vec<AlignmentModelBase>) {
        let restored: Vec<bool> = vec![true; replacement.len()];
        self.model.splice(range.clone(), replacement);
        self.restored.splice(range.start as usize..range.end as usize, restored);
    }
}


#[cfg(test)]
#[path = "../tests/alignment/alignment_model_editor.rs"]
mod tests;