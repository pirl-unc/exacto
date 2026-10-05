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


mod constructor;

use exacto_core::prelude::{BaseQuality, ReadID, ReadPosition};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct AlignmentModel {
    /// Read ID.
    read_id: ReadID,

    /// BAM records.
    records: Vec<AlignmentRecord>,

    /// One entry per base of the read.
    bases: Vec<AlignmentModelBase>,

    /// Events keyed by the pair of read positions they join.
    /// `read_position_1 <= read_position_2`.
    events: HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,

    /// read position -> the other read position it pairs with in `events`.
    events_index: HashMap<ReadPosition, ReadPosition>
}

impl AlignmentModel {
    pub fn get_read_id(&self) -> ReadID {
        self.read_id
    }

    pub fn get_base(&self, read_position: ReadPosition) -> &AlignmentModelBase {
        &self.bases[read_position as usize]
    }

    pub(super) fn get_base_mut(&mut self, read_position: ReadPosition) -> &mut AlignmentModelBase {
        &mut self.bases[read_position as usize]
    }

    pub fn get_bases(&self) -> &[AlignmentModelBase] {
        &self.bases
    }

    pub(super) fn get_bases_mut(&mut self) -> &mut [AlignmentModelBase] {
        &mut self.bases
    }

    pub fn get_events(&self) -> &HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent> {
        &self.events
    }

    pub fn get_event(&self, p1: ReadPosition, p2: ReadPosition) -> Option<&AlignmentModelEvent> {
        self.events.get(&(p1.min(p2), p1.max(p2)))
    }

    pub fn get_event_at(&self, read_position: ReadPosition) -> Option<&AlignmentModelEvent> {
        let other: ReadPosition = *self.events_index.get(&read_position)?;
        self.get_event(read_position, other)
    }

    pub fn get_read_sequence(&self) -> String {
        self.bases.iter().map(|b| b.get_nucleotide().as_str()).collect()
    }

    pub fn get_records(&self) -> &Vec<AlignmentRecord> {
        &self.records
    }

    pub fn iter_base_quality_scores(
        &self,
    ) -> impl ExactSizeIterator<Item = BaseQuality> + DoubleEndedIterator + '_ {
        self.bases
            .iter()
            .map(AlignmentModelBase::get_base_quality)
    }

    /// Terminal soft clips in read coordinates (end exclusive). Internal clips are excluded;
    /// an entirely clipped read has no anchor and therefore no repairable terminal range.
    pub fn terminal_soft_clip_ranges(&self) -> Vec<std::ops::Range<ReadPosition>> {
        let is_clip = |base: &&AlignmentModelBase| *base.get_kind() == AlignmentModelBaseKind::Softclip;
        let first: ReadPosition = self.bases.iter().take_while(is_clip).count() as u32;
        if first == self.num_bases() {
            return Vec::new();
        }
        let last: ReadPosition = self.num_bases() - self.bases.iter().rev().take_while(is_clip).count() as u32;
        [0..first, last..self.num_bases()].into_iter().filter(|range| !range.is_empty()).collect()
    }

    pub fn is_spliced(&self) -> bool {
        self.events.values().any(|event| *event.get_kind() == AlignmentModelEventKind::Splicing)
    }

    pub fn num_bases(&self) -> u32 {
        self.bases.len() as u32
    }

    pub fn num_events(&self) -> u32 {
        self.events.len() as u32
    }

    /// Keep the events `keep` admits; the index follows.
    pub(super) fn retain_events(&mut self, mut keep: impl FnMut(&AlignmentModelEvent) -> bool) {
        self.events.retain(|_, event| keep(event));
        self.reindex_events();
    }

    /// Replace a read range with `replacement` and re-establish everything keyed by read
    /// position. An event with an end inside the range is dropped; a record is clipped to the
    /// bases that survive and dropped when none do.
    pub(super) fn splice(&mut self, range: std::ops::Range<ReadPosition>, replacement: Vec<AlignmentModelBase>) {
        assert!(range.start <= range.end && range.end <= self.num_bases());
        if range.is_empty() && replacement.is_empty() {
            return;
        }
        let delta: i64 = replacement.len() as i64 - (range.end - range.start) as i64;
        let remap = |position: ReadPosition| -> Option<ReadPosition> {
            if position < range.start {
                Some(position)
            } else if position >= range.end {
                Some((position as i64 + delta) as ReadPosition)
            } else {
                None
            }
        };

        let events: Vec<((ReadPosition, ReadPosition), AlignmentModelEvent)> = self.events
            .drain()
            .filter_map(|(_, mut event)| {
                let previous: ReadPosition = remap(event.get_prev_read_position())?;
                let next: ReadPosition = remap(event.get_next_read_position())?;
                event.set_read_positions(previous, next);
                Some(((previous, next), event))
            })
            .collect();
        self.events = events.into_iter().collect();
        self.reindex_events();

        self.records.retain_mut(|record| {
            let first: ReadPosition = if range.contains(&record.read_start) { range.end } else { record.read_start };
            let last: ReadPosition = if range.contains(&record.read_end) {
                let Some(last) = range.start.checked_sub(1) else { return false; };
                last
            } else {
                record.read_end
            };
            if first > last {
                return false;
            }
            record.read_start = remap(first).unwrap();
            record.read_end = remap(last).unwrap();
            true
        });

        self.bases.splice(range.start as usize..range.end as usize, replacement);
        for (position, base) in self.bases.iter_mut().enumerate() {
            base.set_read_position(position as ReadPosition);
        }
    }

    /// Rebuild `events_index` from `events` in key order, so a read position shared by two
    /// events resolves to the later one.
    fn reindex_events(&mut self) {
        let mut keys: Vec<(ReadPosition, ReadPosition)> = self.events.keys().copied().collect();
        keys.sort_unstable();
        self.events_index = keys.into_iter().flat_map(|(a, b)| [(a, b), (b, a)]).collect();
    }
}

impl Clone for AlignmentModel {
    fn clone(&self) -> Self {
        AlignmentModel {
            read_id: self.read_id,
            records: self.records.clone(),
            bases: self.bases.clone(),
            events: self.events.clone(),
            events_index: self.events_index.clone()
        }
    }
}


#[cfg(test)]
#[path = "../../tests/alignment/alignment_model.rs"]
mod tests;