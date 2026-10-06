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
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;
use crate::transcript::transcript_model_base_annotation::TranscriptModelBaseAnnotation;


#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TranscriptModelAnnotation {
    /// Sorted by read start; runs never overlap.
    bases: Vec<AnnotatedReadRun>,

    /// HashMap<(read position 1, read position 2), TranscriptModelEventAnnotation>
    events: HashMap<(ReadPosition, ReadPosition), TranscriptModelEventAnnotation>
}

impl TranscriptModelAnnotation {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_skipped_reference_base(&mut self, p1: ReadPosition, p2: ReadPosition, reference_base: ReferenceBase) {
        self.events
            .entry((p1.min(p2), p1.max(p2)))
            .or_insert_with(TranscriptModelEventAnnotation::new)
            .add_skipped_reference_base(reference_base);
    }

    /// Add a run of skipped reference positions to an event. The run holds no position the
    /// event already holds.
    pub fn add_skipped_reference_run(&mut self, p1: ReadPosition, p2: ReadPosition, run: SkippedReferenceRun) {
        self.events
            .entry((p1.min(p2), p1.max(p2)))
            .or_insert_with(TranscriptModelEventAnnotation::new)
            .add_skipped_reference_run(run);
    }

    /// Index of the run holding `read_position`, if any.
    fn run_index(&self, read_position: ReadPosition) -> Option<usize> {
        let i: usize = self.bases.partition_point(|run| run.read_start <= read_position);
        (i > 0 && read_position <= self.bases[i - 1].read_end).then(|| i - 1)
    }

    pub fn get_base(&self, read_position: ReadPosition) -> Option<&TranscriptModelBaseAnnotation> {
        self.run_index(read_position).map(|i| &self.bases[i].annotation)
    }

    pub fn get_base_context(&self, read_position: ReadPosition) -> Option<&AlignmentModelBaseContext> {
        self.get_base(read_position).map(|a| &a.context)
    }

    pub fn get_event(&self, read_position_1: ReadPosition, read_position_2: ReadPosition) -> Option<&TranscriptModelEventAnnotation> {
        self.events.get(&(read_position_1.min(read_position_2), read_position_1.max(read_position_2)))
    }

    /// Every annotated read position, ascending.
    pub fn get_read_positions(&self) -> impl Iterator<Item = ReadPosition> + '_ {
        self.bases.iter().flat_map(|run| run.read_start..=run.read_end)
    }

    pub fn get_reference_transcript_ids(&self) -> HashSet<ReferenceTranscriptID> {
        self.bases
            .iter()
            .filter_map(|run| run.annotation.reference_transcript_id.clone())
            .collect()
    }

    /// The read positions annotated with `reference_transcript_id`, ascending.
    pub fn get_reference_transcript_read_positions(
        &self,
        reference_transcript_id: &str
    ) -> Vec<ReadPosition> {
        self.bases
            .iter()
            .filter(|run| run.annotation.reference_transcript_id.as_deref() == Some(reference_transcript_id))
            .flat_map(|run| run.read_start..=run.read_end)
            .collect()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty() && self.events.is_empty()
    }

    /// Annotate the read positions `read_start..=read_end`, which lie past every position
    /// annotated so far.
    pub fn push_bases(&mut self, read_start: ReadPosition, read_end: ReadPosition, annotation: TranscriptModelBaseAnnotation) {
        assert!(read_start <= read_end);
        if let Some(last) = self.bases.last_mut() {
            assert!(last.read_end < read_start, "Read positions {}-{} are not past the annotated positions.", read_start, read_end);
            if last.read_end + 1 == read_start && last.annotation == annotation {
                last.read_end = read_end;
                return;
            }
        }
        self.bases.push(AnnotatedReadRun {
            read_start,
            read_end,
            annotation
        });
    }

    pub fn set_base(&mut self, read_position: ReadPosition, annotation: TranscriptModelBaseAnnotation) {
        // Step 1. Carve the position out of the run that holds it, if the annotation changes.
        if let Some(i) = self.run_index(read_position) {
            if self.bases[i].annotation == annotation {
                return;
            }
            let run: AnnotatedReadRun = self.bases.remove(i);
            let mut pieces: Vec<AnnotatedReadRun> = Vec::new();
            if run.read_start < read_position {
                pieces.push(AnnotatedReadRun {
                    read_start: run.read_start,
                    read_end: read_position - 1,
                    annotation: run.annotation.clone()
                });
            }
            if read_position < run.read_end {
                pieces.push(AnnotatedReadRun {
                    read_start: read_position + 1,
                    read_end: run.read_end,
                    annotation: run.annotation
                });
            }
            self.bases.splice(i..i, pieces);
        }

        // Step 2. Insert a one-base run, then coalesce with the run after and the run before.
        let i: usize = self.bases.partition_point(|run| run.read_start < read_position);
        self.bases.insert(i, AnnotatedReadRun {
            read_start: read_position,
            read_end: read_position,
            annotation
        });
        if i + 1 < self.bases.len()
            && self.bases[i].read_end + 1 == self.bases[i + 1].read_start
            && self.bases[i].annotation == self.bases[i + 1].annotation
        {
            let next: AnnotatedReadRun = self.bases.remove(i + 1);
            self.bases[i].read_end = next.read_end;
        }
        if i > 0
            && self.bases[i - 1].read_end + 1 == self.bases[i].read_start
            && self.bases[i - 1].annotation == self.bases[i].annotation
        {
            let current: AnnotatedReadRun = self.bases.remove(i);
            self.bases[i - 1].read_end = current.read_end;
        }
    }

    pub fn set_event(&mut self, p1: ReadPosition, p2: ReadPosition, annotation: TranscriptModelEventAnnotation) {
        self.events.insert((p1.min(p2), p1.max(p2)), annotation);
    }

    pub fn set_event_context(&mut self, p1: ReadPosition, p2: ReadPosition, context: AlignmentModelEventContext) {
        self.events
            .entry((p1.min(p2), p1.max(p2)))
            .or_insert_with(TranscriptModelEventAnnotation::new)
            .context = Some(context);
    }

    pub fn set_boundary_event(&mut self, read_position: ReadPosition) {
        self.set_event_context(
            read_position,
            read_position,
            AlignmentModelEventContext::NonCanonicalSplicing
        );
    }
}


#[derive(Clone, Debug, Serialize, Deserialize)]
struct AnnotatedReadRun {
    read_start: ReadPosition,
    read_end: ReadPosition,
    annotation: TranscriptModelBaseAnnotation
}


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SkippedReferenceRun {
    pub reference_chromosome_id: ReferenceChromosomeID,
    pub reference_start: ReferencePosition,
    pub reference_end: ReferencePosition,
    pub reference_strand: Strand,
    pub reference_gene_id: Option<ReferenceGeneID>,
    pub reference_transcript_id: Option<ReferenceTranscriptID>,
    pub reference_exon_id: Option<ReferenceExonID>,
    pub sequence: Box<str>
}

impl SkippedReferenceRun {
    pub(super) fn from_reference_base(base: &ReferenceBase) -> Self {
        Self {
            reference_chromosome_id: base.reference_chromosome_id,
            reference_start: base.reference_position,
            reference_end: base.reference_position,
            reference_strand: base.reference_strand.clone(),
            reference_gene_id: base.reference_gene_id.clone(),
            reference_transcript_id: base.reference_transcript_id.clone(),
            reference_exon_id: base.reference_exon_id.clone(),
            sequence: base.reference_nucleotide.as_str().into()
        }
    }
    
    pub(super) fn key(&self) -> (ReferenceChromosomeID, ReferencePosition, u8) {
        let strand_rank: u8 = match self.reference_strand {
            Strand::Forward => 0,
            Strand::Reverse => 1,
            Strand::Both => 2,
            Strand::Unknown => 3
        };
        (self.reference_chromosome_id, self.reference_start, strand_rank)
    }

    pub(super) fn contains(&self, chromosome_id: ReferenceChromosomeID, position: ReferencePosition, strand: &Strand) -> bool {
        self.reference_chromosome_id == chromosome_id
            && self.reference_strand == *strand
            && self.reference_start <= position
            && position <= self.reference_end
    }

    pub(super) fn continues_into(&self, next: &Self) -> bool {
        self.reference_chromosome_id == next.reference_chromosome_id
            && self.reference_end + 1 == next.reference_start
            && self.reference_strand == next.reference_strand
            && self.reference_gene_id == next.reference_gene_id
            && self.reference_transcript_id == next.reference_transcript_id
            && self.reference_exon_id == next.reference_exon_id
    }

    pub(super) fn absorb(&mut self, next: Self) {
        self.reference_end = next.reference_end;
        self.sequence = format!("{}{}", self.sequence, next.sequence).into();
    }

    pub fn get_bases(&self) -> Vec<ReferenceBase> {
        self.sequence
            .chars()
            .enumerate()
            .map(|(offset, nucleotide)| ReferenceBase::new(
                self.reference_chromosome_id,
                self.reference_start + offset as u32,
                Nucleotide::from_str(&nucleotide.to_string()).expect("Skipped run holds a non-nucleotide character."),
                self.reference_strand.clone(),
                self.reference_gene_id.clone(),
                self.reference_transcript_id.clone(),
                self.reference_exon_id.clone()
            ))
            .collect()
    }
}


#[cfg(test)]
#[path = "../tests/transcript/transcript_model_annotation.rs"]
mod tests;