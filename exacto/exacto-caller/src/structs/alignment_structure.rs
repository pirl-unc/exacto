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


use bimap::BiMap;
use exacto_core::prelude::*;
use itertools::Itertools;
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct AlignmentStructure {
    read_id: usize,

    /// A vector of read sequence bases.
    bases: Vec<AlignmentStructureBase>,
    
    /// A map between (read_position_1, read_position_2) and its alignment event.
    /// Note that `read_position_1` is smaller than or equal to `read_position_2`.
    events: HashMap<(u32, u32), AlignmentStructureEvent>,

    /// Events index. The key is a read position and the value is another read position
    /// that together appear as a key in `events`.
    events_index: HashMap<u32, u32>
}

/// API methods
impl AlignmentStructure {
    pub fn new(read_id: usize) -> Self {
        Self {
            read_id,
            bases: Vec::new(),
            events: HashMap::new(),
            events_index: HashMap::new()
        }
    }

    pub fn add_base(&mut self, base: AlignmentStructureBase) {
        self.bases.push(base);
    }

    pub fn add_event(&mut self, event: AlignmentStructureEvent) {
        // If there already exists a deletion and a splicing is being added,
        // replace the event with the splicing
        if self.events.contains_key(&(event.get_prev_read_position(), event.get_next_read_position())) {
            let existing_event: &AlignmentStructureEvent = self.get_event(event.get_prev_read_position(), event.get_next_read_position());
            if *event.get_kind() == AlignmentStructureEventKind::Splicing &&
                *existing_event.get_kind() == AlignmentStructureEventKind::Deletion {
                // Make sure the events index already contains the event read positions
                assert!(self.events_index.contains_key(&event.get_prev_read_position()));
                assert!(self.events_index.contains_key(&event.get_next_read_position()));

                // Replace the event
                self.events.insert(
                    (event.get_prev_read_position(), event.get_next_read_position()),
                    event
                );

                return;
            }
        }

        // If there already exists a splicing and a deletion is being added,
        // do not add the event (i.e. keep the splicing)
        if self.events.contains_key(&(event.get_prev_read_position(), event.get_next_read_position())) {
            let existing_event: &AlignmentStructureEvent = self.get_event(event.get_prev_read_position(), event.get_next_read_position());
            if *event.get_kind() == AlignmentStructureEventKind::Deletion &&
                *existing_event.get_kind() == AlignmentStructureEventKind::Splicing {
                // Make sure the events index already contains the event read positions
                assert!(self.events_index.contains_key(&event.get_prev_read_position()));
                assert!(self.events_index.contains_key(&event.get_next_read_position()));

                // Do nothing
                return;
            }
        }

        // Otherwise, first index the event
        self.events_index.insert(event.get_prev_read_position(), event.get_next_read_position());
        self.events_index.insert(event.get_next_read_position(), event.get_prev_read_position());

        // Then, add the event
        self.events.insert(
            (event.get_prev_read_position(), event.get_next_read_position()),
            event
        );
    }

    pub fn contextualize(
        &mut self,
        read_name: &str,
        reference_transcript_sequences: &Vec<&ReferenceTranscriptSequence>,
        gene_annotator: &(impl GeneAnnotator + Sync),
        chromosome_names_map: &BiMap<Box<str>, u16>
    ) {
        // Step 1. Make sure the reference transcript sequences have unique reference gene IDs
        let mut reference_gene_ids: HashSet<Box<str>> = HashSet::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            let inserted: bool = reference_gene_ids.insert(reference_transcript_sequence.get_gene_id().into());
            if !inserted {
                panic!("Only 1 reference transcript per gene ID is allowed.");
            }
        }

        // Step 2. Make sure each reference transcript sequence has at least 1 base that
        // overlaps with one of the self.bases
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            let has_overlap = self.get_bases().iter()
                .filter(|base| matches!(
                    base.get_kind(),
                    AlignmentStructureBaseKind::Match |
                    AlignmentStructureBaseKind::Mismatch |
                    AlignmentStructureBaseKind::Insertion
                ))
                .any(|base| {
                    reference_transcript_sequence.get_bases().iter().any(|base_reference| {
                        base.get_reference_chromosome_id().unwrap() == base_reference.reference_chromosome_id
                            && base.get_reference_position().unwrap() == base_reference.reference_position
                            && base.get_reference_strand().as_ref().unwrap().clone() == base_reference.reference_strand
                    })
                });
            assert!(
                has_overlap,
                "Read name: {}. None of the ReferenceTranscriptSequence bases for {} overlaps with any of the AlignmentStructure bases.",
                read_name, reference_transcript_sequence.get_transcript_id()
            );
        }

        // Step 3. Index the positions of the reference transcript sequences
        let mut reference_transcripts_positions_map: HashMap<(u16, u32, Strand), &ReferenceTranscriptSequence> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            for base in reference_transcript_sequence.get_bases() {
                reference_transcripts_positions_map.insert(
                    (base.reference_chromosome_id,
                        base.reference_position,
                        base.reference_strand.clone()),
                    reference_transcript_sequence
                );
            }
        }

        // Step 4. Identify the context of each AlignmentStructureBase
        for i in 0..self.get_bases_length() {
            let base: &mut AlignmentStructureBase = self.get_base_mut(i);
            if matches!(base.get_kind(),
                AlignmentStructureBaseKind::Match |
                AlignmentStructureBaseKind::Mismatch |
                AlignmentStructureBaseKind::Insertion) {
                let key: (u16, u32, Strand) = (
                    base.get_reference_chromosome_id().unwrap(),
                    base.get_reference_position().unwrap(),
                    base.get_reference_strand().as_ref().unwrap().clone()
                );
                if reference_transcripts_positions_map.contains_key(&key) {
                    // Base is an exonic base
                    let reference_transcript_sequence: &ReferenceTranscriptSequence = reference_transcripts_positions_map.get(&key).unwrap();
                    base.set_context(AlignmentStructureBaseContext::Exonic);
                    base.set_reference_gene_id(reference_transcript_sequence.get_gene_id().into());
                    base.set_reference_transcript_id(reference_transcript_sequence.get_transcript_id());
                    let reference_transcript_id: &str = reference_transcript_sequence.get_transcript_id();
                    let reference_transcript: &Transcript = gene_annotator.get_transcript(reference_transcript_id).unwrap();
                    let base_reference_position: isize = base.get_reference_position().unwrap() as isize;
                    for exon in reference_transcript.exons.values() {
                        let exon_start: isize = exon.start as isize;
                        let exon_end: isize = exon.end as isize;
                        if overlaps(base_reference_position, base_reference_position, exon_start, exon_end) {
                            base.set_reference_exon_id(&*exon.exon_id);
                            break;
                        }
                    }
                    assert_eq!(base.get_reference_exon_id().is_some(), true);
                } else {
                    // Base is an intronic base
                    for reference_transcript_sequence in reference_transcript_sequences.iter() {
                        if base.get_reference_position().unwrap() >= reference_transcript_sequence.get_transcript_start() &&
                            base.get_reference_position().unwrap() <= reference_transcript_sequence.get_transcript_end() {
                            base.set_context(AlignmentStructureBaseContext::Intronic);
                            base.set_reference_gene_id(reference_transcript_sequence.get_gene_id().into());
                            base.set_reference_transcript_id(reference_transcript_sequence.get_transcript_id());
                            break;
                        }
                    }
                    if base.get_context().is_none() {
                        base.set_context(AlignmentStructureBaseContext::Intergenic);
                    }
                }
            }
        }

        // Step 5. Identify reference transcript introns
        let mut reference_transcripts_introns_map: HashMap<u16, HashSet<(u32, u32)>> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            let introns: Vec<(u16, u32, u32)> = reference_transcript_sequence.get_introns();
            for (chromosome_id, start, end) in introns.iter() {
                reference_transcripts_introns_map
                    .entry(*chromosome_id)
                    .or_insert_with(HashSet::new)
                    .insert((*start, *end));
            }
        }

        // Step 6. Identify the context of each AlignmentStructureEvent
        let event_keys: Vec<(u32, u32)> = self.get_events().keys().cloned().collect();
        for (read_position_1, read_position_2) in event_keys {
            let base_1: AlignmentStructureBase = self.get_base(read_position_1).clone();
            let base_2: AlignmentStructureBase = self.get_base(read_position_2).clone();
            let alignment_event: &AlignmentStructureEvent = self.get_event(read_position_1, read_position_2);
            let chromosome_name_1: Box<str> = chromosome_names_map.get_by_right(&base_1.get_reference_chromosome_id().unwrap()).unwrap().clone();
            let chromosome_name_2: Box<str> = chromosome_names_map.get_by_right(&base_2.get_reference_chromosome_id().unwrap()).unwrap().clone();
            let gene_ids_1: HashSet<Box<str>> = gene_annotator.get_gene_ids_at_locus(&*chromosome_name_1, base_1.get_reference_position().unwrap()).into_iter().collect();
            let gene_ids_2: HashSet<Box<str>> = gene_annotator.get_gene_ids_at_locus(&*chromosome_name_2, base_2.get_reference_position().unwrap()).into_iter().collect();
            let gene_ids_disjoint_count_1: usize = gene_ids_1.difference(&gene_ids_2).count();
            let gene_ids_disjoint_count_2: usize = gene_ids_2.difference(&gene_ids_1).count();
            match *alignment_event.get_kind() {
                AlignmentStructureEventKind::Splicing => {
                    let reference_chromosome_id: u16 = base_1.get_reference_chromosome_id().unwrap();
                    let reference_strand: Strand = base_1.get_reference_strand().as_ref().unwrap().clone();
                    let mut reference_start: u32 = base_1.get_reference_position().unwrap();
                    let mut reference_end: u32 = base_2.get_reference_position().unwrap();
                    if reference_strand == Strand::Forward {
                        reference_start = reference_start + 1;
                        reference_end = reference_end - 1;
                    } else {
                        reference_start = reference_start - 1;
                        reference_end = reference_end + 1;
                    }

                    if reference_transcripts_introns_map.contains_key(&reference_chromosome_id) {
                        let reference_transcripts_introns = reference_transcripts_introns_map.get(&reference_chromosome_id).unwrap();
                        if reference_transcripts_introns.contains(&(reference_start, reference_end)) ||
                            reference_transcripts_introns.contains(&(reference_end, reference_start)) {
                            self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::CanonicalSplicing);
                        } else {
                            if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                                self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                            } else {
                                self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                            }
                        }
                    } else {
                        self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                    }
                },
                AlignmentStructureEventKind::Breakpoint => {
                    if base_1.get_reference_chromosome_id().unwrap() == base_2.get_reference_chromosome_id().unwrap() {
                        let left_bases_set: HashSet<(u16, u32)> = (0..=read_position_1)
                            .map(|i| {
                                let base = self.get_base(i);
                                (base.get_reference_chromosome_id().unwrap(), base.get_reference_position().unwrap())
                            })
                            .collect();
                        let right_bases_set: HashSet<(u16, u32)> = (read_position_2..self.get_bases_length())
                            .map(|i| {
                                let base = self.get_base(i);
                                (base.get_reference_chromosome_id().unwrap(), base.get_reference_position().unwrap())
                            })
                            .collect();
                        if left_bases_set.is_disjoint(&right_bases_set) == false {
                            self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::BackSplicing);
                        } else {
                            if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                                self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                            } else {
                                self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                            }
                        }
                    } else {
                        if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                            self.get_event_mut(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                        }
                    }
                },
                _ => {
                    // Do nothing
                }
            }
        }

        // Step 7. Identify skipped reference transcript bases
        let mut reference_transcripts_bases: HashMap<(u16, u32, &Strand), (&str, &ReferenceBase)> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            for base in reference_transcript_sequence.get_bases() {
                reference_transcripts_bases.insert(
                    (base.reference_chromosome_id, base.reference_position, &base.reference_strand),
                    (reference_transcript_sequence.get_transcript_id(), base)
                );
            }
        }
        let mut alignment_structure_bases: HashSet<(u16, u32, &Strand)> = HashSet::new();
        for base in self.get_bases() {
            if matches!(base.get_kind(),
                AlignmentStructureBaseKind::Match |
                AlignmentStructureBaseKind::Mismatch |
                AlignmentStructureBaseKind::Insertion) {
                alignment_structure_bases.insert((
                    base.get_reference_chromosome_id().unwrap(),
                    base.get_reference_position().unwrap(),
                    base.get_reference_strand().as_ref().unwrap()
                ));
            }
        }
        let mut reference_transcript_bases_skipped: HashSet<(u16, u32, &Strand)> = reference_transcripts_bases
            .keys()
            .cloned()
            .filter(|pos| !alignment_structure_bases.contains(pos))
            .sorted_by_key(|&(_, pos, _)| pos)
            .collect();

        // Step 8. Get the reference positions corresponding to exonic boundaries in the
        // alignment structure bases and sort them
        let mut exon_boundary_positions_map: HashMap<Box<str>, Vec<(u32, u32)>> = HashMap::new();
        for ((read_position_1, read_position_2), event) in self.get_events() {
            let base_1: &AlignmentStructureBase = self.get_base(*read_position_1);
            let base_2: &AlignmentStructureBase = self.get_base(*read_position_2);
            if base_1.get_reference_transcript_id().is_some() {
                let reference_transcript_id: Box<str> = base_1.get_reference_transcript_id().as_ref().unwrap().clone();
                exon_boundary_positions_map
                    .entry(reference_transcript_id)
                    .or_insert(Vec::new())
                    .push((*read_position_1, base_1.get_reference_position().unwrap()));
            }
            if base_2.get_reference_transcript_id().is_some() {
                let reference_transcript_id: Box<str> = base_2.get_reference_transcript_id().as_ref().unwrap().clone();
                exon_boundary_positions_map
                    .entry(reference_transcript_id)
                    .or_insert(Vec::new())
                    .push((*read_position_2, base_2.get_reference_position().unwrap()));
            }
        }
        for reference_transcript_id in self.get_reference_transcript_ids() {
            // Include the first and last base of the reference transcript
            // Make sure they are not embedded insertion bases or soft clipped bases
            let mut reference_transcript_bases: Vec<&AlignmentStructureBase> = self.get_reference_transcript_bases(&*reference_transcript_id);
            reference_transcript_bases.sort_by_key(|base| base.get_reference_position().unwrap());
            for base in reference_transcript_bases.iter() {
                if matches!(base.get_kind(),
                    AlignmentStructureBaseKind::Match |
                    AlignmentStructureBaseKind::Mismatch |
                    AlignmentStructureBaseKind::Insertion) {
                    exon_boundary_positions_map
                        .entry(reference_transcript_id.clone())
                        .or_insert(Vec::new())
                        .push(
                            (base.get_read_position(), base.get_reference_position().unwrap())
                        );
                    break;
                }
            }
            for base in reference_transcript_bases.iter().rev() {
                if matches!(base.get_kind(),
                    AlignmentStructureBaseKind::Match |
                    AlignmentStructureBaseKind::Mismatch |
                    AlignmentStructureBaseKind::Insertion) {
                    exon_boundary_positions_map
                        .entry(reference_transcript_id.clone())
                        .or_insert(Vec::new())
                        .push(
                            (base.get_read_position(), base.get_reference_position().unwrap())
                        );
                    break;
                }
            }
        }
        for vec in exon_boundary_positions_map.values_mut() {
            vec.sort_by_key(|&(_, reference_position)| reference_position);
        }

        // Step 9. Identify the closet event for each skipped reference transcript base
        // Map of (chromosome ID, reference position, reference strand) and (read position 1, read position 2)
        let mut events_map: HashMap<(u16, u32, &Strand), (u32, u32)> = HashMap::new();
        for (reference_chromosome_id, reference_position, reference_strand) in reference_transcript_bases_skipped.iter() {
            let reference_transcript_id: Box<str> = reference_transcripts_bases
                .get(&(*reference_chromosome_id, *reference_position, reference_strand))
                .unwrap()
                .0
                .into();

            // Identify the closest base
            let vec: &Vec<(u32, u32)> = exon_boundary_positions_map
                .get(&reference_transcript_id)
                .expect(
                    &format!("Missing reference transcript bases for read name {} and reference transcript ID {}", read_name, reference_transcript_id)
                );
            let closest_read_position: u32 = match vec.binary_search_by_key(reference_position, |&(_, base_reference_position)| base_reference_position) {
                Ok(idx) => {
                    vec[idx].0
                },
                Err(idx) => {
                    if idx == 0 {
                        vec[0].0
                    } else if idx == vec.len() {
                        vec[vec.len() - 1].0
                    } else {
                        let read_position_1: u32 = vec[idx-1].0;
                        let read_position_2: u32 = vec[idx].0;
                        let reference_position_1: u32 = vec[idx-1].1;
                        let reference_position_2: u32 = vec[idx].1;
                        if reference_position.abs_diff(reference_position_1) <= reference_position.abs_diff(reference_position_2) {
                            read_position_1
                        } else {
                            read_position_2
                        }
                    }
                }
            };

            // Identify the closest event
            if closest_read_position == 0 || closest_read_position == self.get_bases_length() - 1 {
                // The skipped reference bases lie beyond the boundaries of the alignment structure for the given reference transcript
                let mut event: AlignmentStructureEvent = AlignmentStructureEvent::new(
                    AlignmentStructureEventKind::Boundary,
                    closest_read_position,
                    closest_read_position,
                    GraphOperationType::Mark,
                    GraphOperationType::Mark
                );
                event.set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                self.add_event(event);
                events_map.insert(
                    (*reference_chromosome_id, *reference_position, reference_strand), (closest_read_position, closest_read_position)
                );
            } else if self.events_index.contains_key(&closest_read_position) {
                let closest_read_position_2: u32 = *self.events_index.get(&closest_read_position).unwrap();
                assert!(self.has_event_between(closest_read_position, closest_read_position_2), "Event does not exist between bases {} and {}", closest_read_position, closest_read_position_2);
                events_map.insert(
                    (*reference_chromosome_id, *reference_position, reference_strand), (closest_read_position, closest_read_position_2)
                );
            } else {
                // No event exists at this position (e.g. no splicing in this read).
                // Create a boundary event to anchor the skipped reference base.
                let mut event: AlignmentStructureEvent = AlignmentStructureEvent::new(
                    AlignmentStructureEventKind::Boundary,
                    closest_read_position,
                    closest_read_position,
                    GraphOperationType::Mark,
                    GraphOperationType::Mark
                );
                event.set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                self.add_event(event);
                events_map.insert(
                    (*reference_chromosome_id, *reference_position, reference_strand), (closest_read_position, closest_read_position)
                );
            }
        }

        assert_eq!(events_map.keys().len(), reference_transcript_bases_skipped.len(), "All skipped reference transcripts' bases should be mapped to an event");

        // Step 10. Insert skipped reference transcript bases into the alignment structure event
        for ((chromosome_id, reference_position, reference_strand), (read_position_1, read_position_2)) in events_map.iter() {
            let reference_base: &ReferenceBase = reference_transcripts_bases.get(&(*chromosome_id, *reference_position, reference_strand)).unwrap().1;
            assert!(self.has_event_between(*read_position_1, *read_position_2) == true);
            self.get_event_mut(*read_position_1, *read_position_2).add_skipped_reference_base(reference_base.clone());
        }
    }

    pub fn get_base(&self, read_position: u32) -> &AlignmentStructureBase {
        match self.bases.get(read_position as usize) {
            Some(base) => {
                base
            },
            None => {
                panic!(
                    "Invalid read_position: {}\n
                     Read length: {}\n
                     First base: {:?}",
                    read_position,
                    self.bases.len(),
                    self.bases.first().unwrap()
                );
            }
        }
    }

    pub fn get_bases(&self) -> &Vec<AlignmentStructureBase> {
        &self.bases
    }

    pub fn get_bases_length(&self) -> u32 {
        self.bases.len() as u32
    }
    
    pub fn get_event(&self,
        read_position_1: u32,
        read_position_2: u32
    ) -> &AlignmentStructureEvent {
        if read_position_1 < read_position_2 {
            self.events
                .get(&(read_position_1, read_position_2))
                .expect(&format!("Missing event: {:?}", self))
        } else {
            self.events
                .get(&(read_position_2, read_position_1))
                .expect(&format!("Missing event: {:?}", self))
        }
    }

    pub fn get_events(&self) -> &HashMap<(u32, u32), AlignmentStructureEvent> {
        &self.events
    }

    pub fn get_event_at_read_position(&self, read_position: u32) -> &AlignmentStructureEvent {
        let read_position_2 = self.events_index.get(&read_position).unwrap();
        self.get_event(read_position, *read_position_2)
    }
    
    pub fn get_events_of_kind(&self, kind: AlignmentStructureEventKind) -> Vec<(u32, u32, &AlignmentStructureEvent)> {
        let mut events: Vec<(u32, u32, &AlignmentStructureEvent)> = Vec::new();
        for ((read_position_1, read_position_2), event) in self.get_events() {
            if event.get_kind() == &kind {
                events.push((*read_position_1, *read_position_2, event));
            }
        }
        events
    }

    pub fn get_base_mut(&mut self, read_position: u32) -> &mut AlignmentStructureBase {
        self.bases.get_mut(read_position as usize).unwrap()
    }

    pub fn get_event_mut(
        &mut self,
        read_position_1: u32,
        read_position_2: u32
    ) -> &mut AlignmentStructureEvent {
        if read_position_1 < read_position_2 {
            assert!(
                self.events.contains_key(&(read_position_1, read_position_2)),
                "Event does not exist between bases {} and {}", read_position_1, read_position_2
            );
            self.events.get_mut(&(read_position_1, read_position_2)).unwrap()
        } else {
            assert!(
                self.events.contains_key(&(read_position_2, read_position_1)),
                "Event does not exist between bases {} and {}", read_position_2, read_position_1
            );
            self.events.get_mut(&(read_position_2, read_position_1)).unwrap()
        }
    }
    
    pub fn get_read_id(&self) -> usize {
        self.read_id
    }
    
    pub fn get_read_sequence(&self) -> Box<str> {
        let mut s: String = String::new();
        for alignment_base in self.bases.iter() {
            s.push_str(alignment_base.get_nucleotide().as_str());
        }
        s.into()
    }

    pub fn get_reference_transcript_bases(&self, reference_transcript_id: &str) -> Vec<&AlignmentStructureBase> {
        let mut reference_transcript_bases: Vec<&AlignmentStructureBase> = Vec::new();
        for base in self.bases.iter() {
            if let Some(reference_transcript_id_) = base.get_reference_transcript_id() {
                if &**reference_transcript_id_ == reference_transcript_id {
                    reference_transcript_bases.push(base);
                }
            }
        }
        reference_transcript_bases.sort_by_key(|base| base.get_reference_position().unwrap());
        reference_transcript_bases
    }

    pub fn get_reference_transcript_ids(&self) -> HashSet<Box<str>> {
        let mut reference_transcript_ids: HashSet<Box<str>> = HashSet::new();
        for base in self.bases.iter() {
            if let Some(reference_transcript_id) = base.get_reference_transcript_id() {
                reference_transcript_ids.insert(reference_transcript_id.clone());
            }
        }
        reference_transcript_ids
    }

    pub fn has_event(&self, read_position: u32) -> bool {
        if self.events_index.contains_key(&read_position) {
            true
        } else {
            false
        }
    }

    pub fn has_event_between(
        &self,
        read_position_1: u32,
        read_position_2: u32
    ) -> bool {
        if read_position_1 < read_position_2 {
            self.events.contains_key(&(read_position_1, read_position_2))
        } else {
            self.events.contains_key(&(read_position_2, read_position_1))
        }
    }

    pub fn identify_exons(&self, read_name: &str) -> Vec<TranscriptModelExon> {
        // Step 1. Cluster adjacent bases by reference position
        let is_exonic = |base: &AlignmentStructureBase| matches!(
            base.get_kind(),
            AlignmentStructureBaseKind::Match |
            AlignmentStructureBaseKind::Mismatch |
            AlignmentStructureBaseKind::Insertion
        );
        let mut uf: UnionFind = UnionFind::new();
        for i in 0..self.get_bases_length() {
            let curr_base: &AlignmentStructureBase = self.get_base(i);
            if is_exonic(curr_base) == false {
                continue;
            }
            uf.union(i, i);

            if i == 0 {
                continue;
            }

            let prev_base: &AlignmentStructureBase = self.get_base(i - 1);
            if is_exonic(prev_base) == false {
                continue;
            }

            let should_union: bool = if self.has_event_between(i - 1, i) {
                let event = self.get_event(i - 1, i);
                event.get_kind() == &AlignmentStructureEventKind::Deletion
                    && prev_base.get_reference_chromosome_id() == curr_base.get_reference_chromosome_id()
                    && prev_base.get_reference_strand() == curr_base.get_reference_strand()
            } else {
                matches!((prev_base.get_reference_chromosome_id(), curr_base.get_reference_chromosome_id()), (Some(p), Some(c)) if p == c) &&
                    prev_base.get_reference_strand() == curr_base.get_reference_strand() &&
                    prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1
            };

            if should_union {
                uf.union(i - 1, i);
            }
        }

        // Step 2. Identify exonic boundaries
        let mut exons: Vec<TranscriptModelExon> = Vec::new();
        let mut clusters: Vec<HashSet<u32>> = uf.get_clusters();
        for cluster in clusters.iter() {
            // Sort the read positions
            let mut read_positions: Vec<u32> = cluster.iter().map(|&pos| pos).collect();
            read_positions.sort();

            let mut bases: Vec<&AlignmentStructureBase> = Vec::new();
            for read_position in read_positions.iter() {
                bases.push(self.get_base(*read_position));
            }

            if let Some(reference_chromosome_id) = bases.first().unwrap().get_reference_chromosome_id() {
                let reference_strand: Strand = bases.first().unwrap().get_reference_strand().as_ref().unwrap().clone();
                let reference_start: u32 = bases.iter().map(|base| base.get_reference_position().unwrap()).min().unwrap();
                let reference_end: u32 = bases.iter().map(|base| base.get_reference_position().unwrap()).max().unwrap();
                let read_start_position: u32 = bases.first().unwrap().get_read_position();
                let read_end_position: u32 = bases.last().unwrap().get_read_position();

                // Use 0 as a placeholder — exon_number assigned after sorting
                let exon: TranscriptModelExon = TranscriptModelExon::new(
                    *reference_chromosome_id,
                    reference_start,
                    reference_end,
                    reference_strand,
                    0,
                    read_start_position,
                    read_end_position
                );

                exons.push(exon);
            }
        }

        // Sort by read order, then assign exon numbers
        exons.sort_by_key(|e| e.read_start_position);
        for (i, exon) in exons.iter_mut().enumerate() {
            exon.exon_number = (i + 1) as u16;
        }

        exons
    }

    pub fn identify_introns(
        &self,
        chromosome_names_map: &BiMap<Box<str>, u16>,
        reference_genome_fasta_file: &str
    ) -> Vec<TranscriptModelIntron> {
        let mut introns: Vec<TranscriptModelIntron> = Vec::new();
        let mut intron_number: u16 = 1;
        for ((read_position_1, read_position_2), event) in self.get_events() {
            match event.get_kind() {
                AlignmentStructureEventKind::Splicing => {
                    let base_1: &AlignmentStructureBase = self.get_base(*read_position_1);
                    let base_2: &AlignmentStructureBase = self.get_base(*read_position_2);
                    let reference_chromosome_id: u16 = base_1.get_reference_chromosome_id().unwrap();
                    let reference_chromosome_name: Box<str> = chromosome_names_map.get_by_right(&reference_chromosome_id).unwrap().clone();
                    let reference_strand: Strand = base_1.get_reference_strand().as_ref().unwrap().clone();
                    let (reference_start_position, reference_end_position) = if reference_strand == Strand::Forward {
                        (base_1.get_reference_position().unwrap() + 1, base_2.get_reference_position().unwrap() - 1)
                    } else {
                        (base_2.get_reference_position().unwrap() + 1, base_1.get_reference_position().unwrap() - 1)
                    };

                    // Make sure the bases are on the same reference chromosome
                    assert!(reference_chromosome_id == base_2.get_reference_chromosome_id().unwrap());

                    // Make sure the bases are on the same reference strand
                    assert!(reference_strand == base_2.get_reference_strand().as_ref().unwrap().clone());

                    let sequence_1: Box<str> = get_fasta_sequence(
                        &*reference_chromosome_name,
                        reference_start_position + 1,
                        reference_start_position + 2,
                        reference_genome_fasta_file
                    );

                    let sequence_2: Box<str> = get_fasta_sequence(
                        &*reference_chromosome_name,
                        reference_end_position - 2,
                        reference_end_position - 1,
                        reference_genome_fasta_file
                    );

                    let (donor_splice_site_signal, acceptor_splice_site_signal) = if reference_strand == Strand::Forward {
                        (sequence_1, sequence_2)
                    } else {
                        (reverse_complement(&*sequence_2).into(), reverse_complement(&*sequence_1).into())
                    };

                    let intron: TranscriptModelIntron = TranscriptModelIntron::new(
                        reference_chromosome_id,
                        reference_start_position,
                        reference_end_position,
                        reference_strand,
                        intron_number,
                        &*donor_splice_site_signal,
                        &*acceptor_splice_site_signal,
                        *read_position_1,
                        *read_position_2
                    );

                    intron_number += 1;
                    introns.push(intron);
                },
                _ => {
                    // Do nothing
                }
            }
        }

        introns
    }

    pub fn identify_records(&self) -> Vec<AlignmentStructureRecord> {
        // Step 1. Cluster Match, Mismatch, Insertion bases
        let num_bases: u32 = self.get_bases_length();
        let mut uf_bases: UnionFind = UnionFind::new();
        for i in 0..num_bases {
            if matches!(self.get_base(i).get_kind(),
                AlignmentStructureBaseKind::Match |
                AlignmentStructureBaseKind::Mismatch |
                AlignmentStructureBaseKind::Insertion) {
                uf_bases.union(i, i);
            }
        }
        for i in 0..num_bases {
            if i > 0 {
                let prev_base: &AlignmentStructureBase = self.get_base(i - 1);
                let curr_base: &AlignmentStructureBase = self.get_base(i);
                let is_prev_base_aligned: bool = matches!(prev_base.get_kind(),
                    AlignmentStructureBaseKind::Match |
                    AlignmentStructureBaseKind::Mismatch |
                    AlignmentStructureBaseKind::Insertion
                );
                let is_curr_base_aligned: bool = matches!(curr_base.get_kind(),
                    AlignmentStructureBaseKind::Match |
                    AlignmentStructureBaseKind::Mismatch |
                    AlignmentStructureBaseKind::Insertion
                );
                if is_prev_base_aligned && is_curr_base_aligned {
                    if prev_base.get_context().is_some() && curr_base.get_context().is_some() {
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            *prev_base.get_context().as_ref().unwrap() == *curr_base.get_context().as_ref().unwrap() &&
                            prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                            prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1 &&
                            prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() {
                            uf_bases.union(i - 1, i);
                        }
                    } else {
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                            prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1 &&
                            prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() {
                            uf_bases.union(i - 1, i);
                        }
                    }
                }
            }
        }

        // Step 2. Record bases
        let mut records: Vec<AlignmentStructureRecord> = Vec::new();
        for cluster in uf_bases.get_clusters() {
            let mut read_positions: Vec<u32> = cluster.into_iter().collect();
            read_positions.sort();
            let bases: Vec<&AlignmentStructureBase> = read_positions.iter().map(|&i| self.get_base(i)).collect();
            let first_base: &AlignmentStructureBase = bases.first().unwrap();
            let last_base: &AlignmentStructureBase = bases.last().unwrap();
            let mut sequence: String = String::new();
            let mut base_quality_scores: Vec<u8> = Vec::new();
            for base in bases.iter() {
                assert_eq!(
                    matches!(base.get_kind(),
                        AlignmentStructureBaseKind::Match |
                        AlignmentStructureBaseKind::Mismatch |
                        AlignmentStructureBaseKind::Insertion
                    ),
                    true
                );
                sequence.push_str(base.get_nucleotide().as_str());
                base_quality_scores.push(base.get_base_quality());
            }
            let (base_1, base_2) = if *first_base.get_reference_strand().as_ref().unwrap() == Strand::Forward {
                (first_base, last_base)
            } else {
                (last_base, first_base)
            };
            match first_base.get_kind() {
                AlignmentStructureBaseKind::Match => {
                    let record: AlignmentStructureRecord = AlignmentStructureRecord::new(
                        first_base.get_read_position(),
                        last_base.get_read_position(),
                        sequence.as_str(),
                        base_quality_scores,
                        AlignmentStructureRecordType::Base,
                        AlignmentStructureKind::Base(first_base.get_kind().clone()),
                        first_base.get_context().as_ref().cloned().map(AlignmentStructureContext::Base),
                        base_1.get_reference_chromosome_id().unwrap(),
                        base_1.get_reference_position().unwrap(),
                        GraphOperationType::Include,
                        base_1.get_reference_strand().as_ref().unwrap().clone(),
                        base_1.get_mapping_quality().unwrap_or(0),
                        base_2.get_reference_chromosome_id().unwrap(),
                        base_2.get_reference_position().unwrap(),
                        GraphOperationType::Include,
                        base_2.get_reference_strand().as_ref().unwrap().clone(),
                        base_2.get_mapping_quality().unwrap_or(0),
                        base_1.get_reference_gene_id().clone(),
                        base_1.get_reference_transcript_id().clone(),
                        base_1.get_reference_exon_id().clone(),
                        base_2.get_reference_gene_id().clone(),
                        base_2.get_reference_transcript_id().clone(),
                        base_2.get_reference_exon_id().clone(),
                        None
                    );
                    records.push(record);
                },
                AlignmentStructureBaseKind::Mismatch => {
                    let record: AlignmentStructureRecord = AlignmentStructureRecord::new(
                        first_base.get_read_position(),
                        last_base.get_read_position(),
                        sequence.as_str(),
                        base_quality_scores,
                        AlignmentStructureRecordType::Base,
                        AlignmentStructureKind::Base(first_base.get_kind().clone()),
                        first_base.get_context().as_ref().cloned().map(AlignmentStructureContext::Base),
                        base_1.get_reference_chromosome_id().unwrap(),
                        base_1.get_reference_position().unwrap() - 1,
                        GraphOperationType::Downstream,
                        base_1.get_reference_strand().as_ref().unwrap().clone(),
                        base_1.get_mapping_quality().unwrap_or(0),
                        base_2.get_reference_chromosome_id().unwrap(),
                        base_2.get_reference_position().unwrap() + 1,
                        GraphOperationType::Upstream,
                        base_2.get_reference_strand().as_ref().unwrap().clone(),
                        base_2.get_mapping_quality().unwrap_or(0),
                        base_1.get_reference_gene_id().clone(),
                        base_1.get_reference_transcript_id().clone(),
                        base_1.get_reference_exon_id().clone(),
                        base_2.get_reference_gene_id().clone(),
                        base_2.get_reference_transcript_id().clone(),
                        base_2.get_reference_exon_id().clone(),
                        None
                    );
                    records.push(record);
                },
                AlignmentStructureBaseKind::Insertion => {
                    let record: AlignmentStructureRecord = AlignmentStructureRecord::new(
                        first_base.get_read_position(),
                        last_base.get_read_position(),
                        sequence.as_str(),
                        base_quality_scores,
                        AlignmentStructureRecordType::Base,
                        AlignmentStructureKind::Base(first_base.get_kind().clone()),
                        first_base.get_context().as_ref().cloned().map(AlignmentStructureContext::Base),
                        base_1.get_reference_chromosome_id().unwrap(),
                        base_1.get_reference_position().unwrap(),
                        GraphOperationType::Downstream,
                        base_1.get_reference_strand().as_ref().unwrap().clone(),
                        base_1.get_mapping_quality().unwrap_or(0),
                        base_2.get_reference_chromosome_id().unwrap(),
                        base_2.get_reference_position().unwrap() + 1,
                        GraphOperationType::Upstream,
                        base_2.get_reference_strand().as_ref().unwrap().clone(),
                        base_2.get_mapping_quality().unwrap_or(0),
                        base_1.get_reference_gene_id().clone(),
                        base_1.get_reference_transcript_id().clone(),
                        base_1.get_reference_exon_id().clone(),
                        base_2.get_reference_gene_id().clone(),
                        base_2.get_reference_transcript_id().clone(),
                        base_2.get_reference_exon_id().clone(),
                        None
                    );
                    records.push(record);
                },
                _ => {
                    // Do nothing
                }
            }
        }

        // Step 3. Record events
        for ((read_position_1, read_position_2), event) in self.events.iter() {
            let prev_base: &AlignmentStructureBase = self.get_base(event.get_prev_read_position());
            let next_base: &AlignmentStructureBase = self.get_base(event.get_next_read_position());

            // Get the sequence between the two read positions
            let mut sequence: String = "".to_string();
            let mut base_quality_scores: Vec<u8> = Vec::new();
            if read_position_1 < read_position_2 {
                for k in read_position_1 + 1..=read_position_2 - 1 {
                    let base: &AlignmentStructureBase = self.get_base(k);
                    sequence.push_str(base.get_nucleotide().as_str());
                    base_quality_scores.push(base.get_base_quality());
                }
            }

            let (base_1, base_2, operation_1, operation_2) = if prev_base.get_reference_chromosome_id().unwrap() == next_base.get_reference_chromosome_id().unwrap() {
                if prev_base.get_reference_position().unwrap() < next_base.get_reference_position().unwrap() {
                    (prev_base, next_base, event.get_prev_graph_operation_type(), event.get_next_graph_operation_type())
                } else {
                    (next_base, prev_base, event.get_next_graph_operation_type(), event.get_prev_graph_operation_type())
                }
            } else {
                (prev_base, next_base, event.get_prev_graph_operation_type(), event.get_next_graph_operation_type())
            };

            let record: AlignmentStructureRecord = AlignmentStructureRecord::new(
                prev_base.get_read_position(),
                next_base.get_read_position(),
                sequence.as_str(),
                base_quality_scores,
                AlignmentStructureRecordType::Event,
                AlignmentStructureKind::Event(event.get_kind().clone()),
                event.get_context().as_ref().cloned().map(AlignmentStructureContext::Event),
                base_1.get_reference_chromosome_id().unwrap(),
                base_1.get_reference_position().unwrap(),
                operation_1.clone(),
                base_1.get_reference_strand().as_ref().unwrap().clone(),
                base_1.get_mapping_quality().unwrap_or(0),
                base_2.get_reference_chromosome_id().unwrap(),
                base_2.get_reference_position().unwrap(),
                operation_2.clone(),
                base_2.get_reference_strand().as_ref().unwrap().clone(),
                base_2.get_mapping_quality().unwrap_or(0),
                base_1.get_reference_gene_id().clone(),
                base_1.get_reference_transcript_id().clone(),
                base_1.get_reference_exon_id().clone(),
                base_2.get_reference_gene_id().clone(),
                base_2.get_reference_transcript_id().clone(),
                base_2.get_reference_exon_id().clone(),
                Some(event.get_skipped_reference_bases_clusters())
            );
            records.push(record);
        }

        // Step 4. Sort the records
        records.sort_by(|a, b| {
            a.get_start().cmp(&b.get_start()).then(a.get_end().cmp(&b.get_end()))
        });

        records
    }

    pub fn identify_variant_records(
        &self,
        min_mapping_quality: u16,
        min_base_quality: u8,
        analyte_type: AnalyteType
    ) -> Vec<VariantRecord> {
        // Step 1. Get records
        let records: Vec<AlignmentStructureRecord> = self.identify_records();

        // Step 2. Identify base variant records
        let mut variant_records: Vec<VariantRecord> = self.identify_base_variant_records(
            &records,
            min_mapping_quality,
            min_base_quality,
            analyte_type.clone()
        );

        // Step 3. Identify terminal-softclip insertion variant records
        variant_records.extend(
            self.identify_terminal_softclip_insertions(
                min_mapping_quality,
                min_base_quality,
                analyte_type.clone()
            )
        );

        // Step 4. Identify event variant records
        variant_records.extend(self.identify_event_variant_records(
            &records,
            min_mapping_quality,
            min_base_quality
        ));

        // Step 5. Sort the variant records
        variant_records.sort_by(|a, b| {
            a.get_chromosome_1()
                .cmp(&b.get_chromosome_1())
                .then(a.get_position_1().cmp(&b.get_position_1()))
        });

        variant_records
    }

    pub fn is_spliced(&self) -> bool {
        for event in self.get_events().values() {
            match event.get_kind() {
                AlignmentStructureEventKind::Splicing => {
                    return true;
                },
                _ => {
                    // Do nothing
                }
            }
        }
        false
    }

    pub fn to_dataframe(&self, chromosome_names_map: &BiMap<Box<str>, u16>) -> DataFrame {
        let records_map = self.to_record(chromosome_names_map);
        DataFrame::new(vec![
            Column::from(Series::new("index".into(), records_map.get("index").unwrap())),
            Column::from(Series::new("read_start".into(), records_map.get("read_start").unwrap())),
            Column::from(Series::new("read_end".into(), records_map.get("read_end").unwrap())),
            Column::from(Series::new("sequence".into(), records_map.get("sequence").unwrap())),
            Column::from(Series::new("type".into(), records_map.get("type").unwrap())),
            Column::from(Series::new("kind".into(), records_map.get("kind").unwrap())),
            Column::from(Series::new("context".into(), records_map.get("context").unwrap())),
            Column::from(Series::new("chromosome_1".into(), records_map.get("chromosome_1").unwrap())),
            Column::from(Series::new("position_1".into(), records_map.get("position_1").unwrap())),
            Column::from(Series::new("operation_1".into(), records_map.get("operation_1").unwrap())),
            Column::from(Series::new("strand_1".into(), records_map.get("strand_1").unwrap())),
            Column::from(Series::new("chromosome_2".into(), records_map.get("chromosome_2").unwrap())),
            Column::from(Series::new("position_2".into(), records_map.get("position_2").unwrap())),
            Column::from(Series::new("operation_2".into(), records_map.get("operation_2").unwrap())),
            Column::from(Series::new("strand_2".into(), records_map.get("strand_2").unwrap())),
            Column::from(Series::new("gene_id_1".into(), records_map.get("gene_id_1").unwrap())),
            Column::from(Series::new("transcript_id_1".into(), records_map.get("transcript_id_1").unwrap())),
            Column::from(Series::new("exon_id_1".into(), records_map.get("exon_id_1").unwrap())),
            Column::from(Series::new("gene_id_2".into(), records_map.get("gene_id_2").unwrap())),
            Column::from(Series::new("transcript_id_2".into(), records_map.get("transcript_id_2").unwrap())),
            Column::from(Series::new("exon_id_2".into(), records_map.get("exon_id_2").unwrap())),
            Column::from(Series::new("skipped".into(), records_map.get("skipped").unwrap()))
        ]).unwrap()
    }

    pub fn to_record(&self, chromosome_names_map: &BiMap<Box<str>, u16>) -> HashMap<Box<str>, Vec<AnyValue>> {
        let mut records_map: HashMap<Box<str>, Vec<AnyValue>> = HashMap::new();
        let mut index: u64 = 0;
        for record in self.identify_records().iter() {
            let chromosome_1: String = chromosome_names_map.get_by_right(&record.get_chromosome_1()).unwrap().to_string();
            let chromosome_2: String = chromosome_names_map.get_by_right(&record.get_chromosome_2()).unwrap().to_string();
            records_map.entry("index".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(index));
            records_map.entry("read_start".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_start() as u64));
            records_map.entry("read_end".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_end() as u64));
            records_map.entry("sequence".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_sequence().to_string().as_str().into()));
            records_map.entry("type".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_record_type().as_str().into()));
            records_map.entry("kind".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_kind().as_str().into()));
            records_map.entry("context".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_context().as_ref().map_or("".into(), |k| k.clone().as_str().into())));
            records_map.entry("chromosome_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(chromosome_1.as_str().into()));
            records_map.entry("position_1".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_position_1() as u64));
            records_map.entry("operation_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_operation_1().as_str().into()));
            records_map.entry("strand_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_strand_1().as_str().into()));
            records_map.entry("chromosome_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(chromosome_2.as_str().into()));
            records_map.entry("position_2".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_position_2() as u64));
            records_map.entry("operation_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_operation_2().as_str().into()));
            records_map.entry("strand_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_strand_2().as_str().into()));
            records_map.entry("gene_id_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_gene_id_1().as_ref().map_or("".into(),|k| k.to_string().into())));
            records_map.entry("transcript_id_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_transcript_id_1().as_ref().map_or("".into(), |k| k.to_string().into())));
            records_map.entry("exon_id_1".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_exon_id_1().as_ref().map_or("".into(), |k| k.to_string().into())));
            records_map.entry("gene_id_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_gene_id_2().as_ref().map_or("".into(), |k| k.to_string().into())));
            records_map.entry("transcript_id_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_transcript_id_2().as_ref().map_or("".into(), |k| k.to_string().into())));
            records_map.entry("exon_id_2".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_exon_id_2().as_ref().map_or("".into(), |k| k.to_string().into())));
            records_map.entry("skipped".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_skipped_string(chromosome_names_map).as_str().into()));
            index += 1;
        }
        records_map
    }
}

/// Helper functions
impl AlignmentStructure {

    /// Build a single `Insertion` `VariantRecord` covering read positions
    /// `[start, end_inclusive]`. Mirrors the position/operation conventions
    /// used by `identify_records` for Insertion-kind base clusters:
    ///   position_1 = base.reference_position, op_1 = Downstream
    ///   position_2 = base.reference_position + 1, op_2 = Upstream
    fn build_softclip_insertion_record(
        &self,
        start: u32,
        end_inclusive: u32,
        min_mapping_quality: u16,
        min_base_quality: u8,
    ) -> Option<VariantRecord> {
        let first: &AlignmentStructureBase = self.get_base(start);

        // All terminal-softclip bases set in alignment.rs::identify_breakpoint_variants
        // share chromosome_id / reference_position / reference_strand / mapping_quality,
        // so reading them off the first base is sufficient.
        let mq: u16 = first.get_mapping_quality().unwrap();
        if mq < min_mapping_quality {
            return None;
        }
        let chromosome_id: u16 = first.get_reference_chromosome_id().unwrap();
        let ref_position: u32 = first.get_reference_position().unwrap();
        let ref_strand: Strand = first.get_reference_strand().as_ref().unwrap().clone();

        // Filter the inserted sequence by base quality.
        let mut sequence: String = String::new();
        for i in start..=end_inclusive {
            let b: &AlignmentStructureBase = self.get_base(i);
            if b.get_base_quality() >= min_base_quality {
                sequence.push_str(b.get_nucleotide().as_str());
            }
        }
        if sequence.is_empty() {
            return None;
        }

        let graph_operation: GraphOperation = GraphOperation::new(
            chromosome_id,
            ref_position,
            ref_strand.clone(),
            GraphOperationType::Downstream,
            chromosome_id,
            ref_position + 1,
            ref_strand,
            GraphOperationType::Upstream,
            sequence.into(),
            VariantType::Insertion,
        );

        Some(VariantRecord::new(
            self.read_id,
            start,
            end_inclusive,
            graph_operation,
        ))
    }

    /// Detect runs of `Softclip`-kind bases that touch the read termini
    /// (start at read position 0 or end at the last read position) and emit
    /// each as a single `Insertion` `VariantRecord`. Mid-read soft-clip runs
    /// (between split alignments) are ignored — those already get a
    /// `Breakpoint` event added by `identify_breakpoint_variants`.
    fn identify_terminal_softclip_insertions(
        &self,
        min_mapping_quality: u16,
        min_base_quality: u8,
        analyte_type: AnalyteType
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        let n: u32 = self.get_bases_length();
        if n == 0 {
            return variant_records;
        }

        // Leading run: contiguous Softclip bases starting at read position 0.
        let mut leading_end: u32 = 0; // exclusive
        while leading_end < n
            && *self.get_base(leading_end).get_kind() == AlignmentStructureBaseKind::Softclip
        {
            leading_end += 1;
        }
        if leading_end > 0 {
            if let Some(vr) = self.build_softclip_insertion_record(
                0,
                leading_end - 1,
                min_mapping_quality,
                min_base_quality,
            ) {
                variant_records.push(vr);
            }
        }

        // Trailing run: contiguous Softclip bases ending at read position n-1.
        let mut trailing_start: u32 = n; // inclusive
        while trailing_start > 0
            && *self.get_base(trailing_start - 1).get_kind() == AlignmentStructureBaseKind::Softclip
        {
            trailing_start -= 1;
        }
        // Skip if there is no trailing run, or if the trailing run overlaps the
        // leading run we already emitted (entire read is one softclip span).
        if trailing_start < n && trailing_start > leading_end {
            if let Some(vr) = self.build_softclip_insertion_record(
                trailing_start,
                n - 1,
                min_mapping_quality,
                min_base_quality,
            ) {
                variant_records.push(vr);
            }
        }

        // Exclude template-switching artifacts if the analyte type is RNA
        if analyte_type == AnalyteType::RNA {
            variant_records.retain(|vr| {
                let sequence: String = vr.get_sequence().to_uppercase();
                if sequence.len() <= 3 {
                    if sequence.contains("CC") ||
                        sequence.contains("GG") ||
                        sequence.contains("GG") ||
                        sequence.contains("CC") ||
                        sequence.contains("GG") {
                        return false;
                    } else {
                        return true;
                    }
                } else {
                    return true;
                }
            });
        }

        variant_records
    }

    /// Identifies base variant records.
    ///
    /// This function identifies the following variant types:
    /// - Single-nucleotide variant
    /// - Multi-nucleotide variant
    /// - Insertion
    /// - Cryptic exon
    /// - Intron retention
    /// - UTR extension
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality.
    /// * `min_base_quality`: Minimum base quality.
    ///
    /// # Returns
    /// * Vector of VariantRecord objects.
    fn identify_base_variant_records(
        &self,
        records: &Vec<AlignmentStructureRecord>,
        min_mapping_quality: u16,
        min_base_quality: u8,
        analyte_type: AnalyteType
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();

        // Step 1. Identify variant records based on the AlignmentStructureBase kinds
        for i in 0..records.len() {
            let curr_record: &AlignmentStructureRecord = records.get(i).unwrap();

            if *curr_record.get_record_type() == AlignmentStructureRecordType::Event {
                continue;
            }

            // Check the mapping quality scores
            if curr_record.get_mapping_quality_1() < min_mapping_quality || curr_record.get_mapping_quality_2() < min_mapping_quality {
                continue;
            }

            match curr_record.get_kind() {
                AlignmentStructureKind::Base(AlignmentStructureBaseKind::Mismatch) => {
                    // Exclude SNVs or MNVs at the first or last base of the read
                    // if the analyte type is RNA
                    if analyte_type == AnalyteType::RNA &&
                        (curr_record.get_start() == 0 || curr_record.get_end() == self.get_bases_length() - 1) {
                        continue;
                    }
                    let sequence: String = curr_record
                        .get_sequence()
                        .chars()
                        .zip(curr_record.get_base_quality_scores())
                        .filter(|(_, &q)| q >= min_base_quality)
                        .map(|(c, _)| c)
                        .collect();

                    if sequence.len() == 0 {
                        continue;
                    }

                    let decrement: u32 = sequence.len().abs_diff(curr_record.get_sequence().len()) as u32;

                    let variant_type = if sequence.len() == 1 {
                        VariantType::SingleNucleotideVariant
                    } else {
                        VariantType::MultiNucleotideVariant
                    };

                    let graph_operation: GraphOperation = GraphOperation::new(
                        curr_record.get_chromosome_1(),
                        curr_record.get_position_1(),
                        curr_record.get_strand_1().clone(),
                        curr_record.get_operation_1().clone(),
                        curr_record.get_chromosome_2(),
                        curr_record.get_position_2() - decrement,
                        curr_record.get_strand_2().clone(),
                        curr_record.get_operation_2().clone(),
                        sequence.into(),
                        variant_type
                    );

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            curr_record.get_start(),
                            curr_record.get_end(),
                            graph_operation
                        )
                    );
                },
                AlignmentStructureKind::Base(AlignmentStructureBaseKind::Insertion) => {
                    let sequence: String = curr_record
                        .get_sequence()
                        .chars()
                        .zip(curr_record.get_base_quality_scores())
                        .filter(|(_, &q)| q >= min_base_quality)
                        .map(|(c, _)| c)
                        .collect();

                    if sequence.len() == 0 {
                        continue;
                    }

                    let graph_operation: GraphOperation = GraphOperation::new(
                        curr_record.get_chromosome_1(),
                        curr_record.get_position_1(),
                        curr_record.get_strand_1().clone(),
                        curr_record.get_operation_1().clone(),
                        curr_record.get_chromosome_2(),
                        curr_record.get_position_2(),
                        curr_record.get_strand_2().clone(),
                        curr_record.get_operation_2().clone(),
                        sequence.into(),
                        VariantType::Insertion
                    );

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            curr_record.get_start(),
                            curr_record.get_end(),
                            graph_operation
                        )
                    );
                },
                _ => {
                    // Do nothing
                }
            }
        }

        // Return the current vector of VariantRecord objects for DNA variant calling
        if analyte_type == AnalyteType::DNA {
            return variant_records;
        }

        // Step 2. Get a set of all included read bases
        let mut included_read_bases: HashSet<(u32, u32)> = HashSet::new();
        for variant_record in variant_records.iter() {
            included_read_bases.insert(
                (variant_record.read_position_1, variant_record.read_position_2)
            );
        }

        // Step 3. Identify variant records based on the AlignmentStructureBase contexts - RNA variant calling
        let mut uf: UnionFind = UnionFind::new();
        for i in 0..records.len() {
            let record: &AlignmentStructureRecord = records.get(i).unwrap();
            if *record.get_record_type() == AlignmentStructureRecordType::Base &&
                *record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) {
                uf.union(i as u32, i as u32);
            }
        }
        for i in 1..records.len() {
            let prev_record: &AlignmentStructureRecord = records.get(i - 1).unwrap();
            let curr_record: &AlignmentStructureRecord = records.get(i).unwrap();
            match (prev_record.get_record_type(), curr_record.get_record_type()) {
                (AlignmentStructureRecordType::Event, AlignmentStructureRecordType::Event) => {
                    continue;
                },
                (AlignmentStructureRecordType::Event, AlignmentStructureRecordType::Base) => {
                    if *curr_record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                        *prev_record.get_kind() == AlignmentStructureKind::Event(AlignmentStructureEventKind::Deletion) &&
                        prev_record.get_chromosome_1() == curr_record.get_chromosome_1() &&
                        prev_record.get_chromosome_2() == curr_record.get_chromosome_2() &&
                        prev_record.get_strand_1() == curr_record.get_strand_1() &&
                        prev_record.get_strand_2() == curr_record.get_strand_2() &&
                        (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1 ||
                         prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                        let base_1: &AlignmentStructureBase = self.get_base(prev_record.get_start());
                        let base_2: &AlignmentStructureBase = self.get_base(prev_record.get_end());
                        if *base_1.get_context().as_ref().unwrap() != AlignmentStructureBaseContext::Exonic &&
                            *base_1.get_context().as_ref().unwrap() == *base_2.get_context().as_ref().unwrap() &&
                            *base_1.get_context().as_ref().unwrap() == *curr_record.get_context().as_ref().unwrap().as_base().unwrap() {
                            uf.union(i as u32 - 1, i as u32);
                        }
                    }
                },
                (AlignmentStructureRecordType::Base, AlignmentStructureRecordType::Event) => {
                    if *prev_record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                        *curr_record.get_kind() == AlignmentStructureKind::Event(AlignmentStructureEventKind::Deletion) &&
                        prev_record.get_chromosome_1() == curr_record.get_chromosome_1() &&
                        prev_record.get_chromosome_2() == curr_record.get_chromosome_2() &&
                        prev_record.get_strand_1() == curr_record.get_strand_1() &&
                        prev_record.get_strand_2() == curr_record.get_strand_2() &&
                        (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1 ||
                         prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                        let base_1: &AlignmentStructureBase = self.get_base(curr_record.get_start());
                        let base_2: &AlignmentStructureBase = self.get_base(curr_record.get_end());
                        if *base_1.get_context().as_ref().unwrap() != AlignmentStructureBaseContext::Exonic &&
                            *base_1.get_context().as_ref().unwrap() == *base_2.get_context().as_ref().unwrap() &&
                            *base_1.get_context().as_ref().unwrap() == *prev_record.get_context().as_ref().unwrap().as_base().unwrap() {
                            uf.union(i as u32 - 1, i as u32);
                        }
                    }
                }
                (AlignmentStructureRecordType::Base, AlignmentStructureRecordType::Base) => {
                    if *prev_record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                        *curr_record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                        *prev_record.get_context().as_ref().unwrap().as_base().unwrap() == *curr_record.get_context().as_ref().unwrap().as_base().unwrap() &&
                        prev_record.get_chromosome_1() == curr_record.get_chromosome_1() &&
                        prev_record.get_chromosome_2() == curr_record.get_chromosome_2() &&
                        prev_record.get_strand_1() == curr_record.get_strand_1() &&
                        prev_record.get_strand_2() == curr_record.get_strand_2() &&
                        (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1 ||
                         prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                        uf.union(i as u32 - 1, i as u32);
                    }
                }
            }
        }
        for cluster in uf.get_clusters().iter() {
            // Sort the record positions
            let mut record_indices: Vec<u32> = cluster.iter().map(|&pos| pos).collect();
            record_indices.sort();

            let first_record: &AlignmentStructureRecord = records.get(*record_indices.first().unwrap() as usize).unwrap();
            let last_record: &AlignmentStructureRecord = records.get(*record_indices.last().unwrap() as usize).unwrap();

            let prev_record: Option<&AlignmentStructureRecord> = if *record_indices.first().unwrap() > 0u32  {
                records.get(record_indices[0] as usize - 1)
            } else {
                None
            };

            let next_record: Option<&AlignmentStructureRecord> = if *record_indices.last().unwrap() < records.len() as u32 - 1  {
                records.get(*record_indices.last().unwrap() as usize + 1)
            } else {
                None
            };

            let mut variant_type: VariantType = VariantType::CrypticExon;
            if prev_record.is_some() {
                let prev_record_: &AlignmentStructureRecord = prev_record.unwrap();
                if *prev_record_.get_record_type() == AlignmentStructureRecordType::Base &&
                    *prev_record_.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                    prev_record_.get_chromosome_1() == first_record.get_chromosome_1() &&
                    prev_record_.get_chromosome_2() == first_record.get_chromosome_2() &&
                    prev_record_.get_strand_1() == first_record.get_strand_1() &&
                    prev_record_.get_strand_2() == first_record.get_strand_2() &&
                    (prev_record_.get_position_1().abs_diff(first_record.get_position_2()) <= 1 ||
                     prev_record_.get_position_2().abs_diff(first_record.get_position_1()) <= 1) {
                    if *first_record.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Intronic) {
                        variant_type = VariantType::IntronRetention;
                    }
                    if *first_record.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Intergenic) {
                        variant_type = VariantType::UTRExtension;
                    }
                }
            }
            if next_record.is_some() {
                let next_record_: &AlignmentStructureRecord = next_record.unwrap();
                if *next_record_.get_record_type() == AlignmentStructureRecordType::Base &&
                    *next_record_.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) &&
                    next_record_.get_chromosome_1() == last_record.get_chromosome_1() &&
                    next_record_.get_chromosome_2() == last_record.get_chromosome_2() &&
                    next_record_.get_strand_1() == last_record.get_strand_1() &&
                    next_record_.get_strand_2() == last_record.get_strand_2() &&
                    (next_record_.get_position_1().abs_diff(last_record.get_position_2()) <= 1 ||
                        next_record_.get_position_2().abs_diff(last_record.get_position_1()) <= 1) {
                    if *last_record.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Intronic) {
                        variant_type = VariantType::IntronRetention;
                    }
                    if *last_record.get_context().as_ref().unwrap() == AlignmentStructureContext::Base(AlignmentStructureBaseContext::Intergenic) {
                        variant_type = VariantType::UTRExtension;
                    }
                }
            }

            assert_eq!(first_record.get_chromosome_1() == last_record.get_chromosome_1(), true, "The first and last records must have the same chromosome 1.");
            assert_eq!(first_record.get_chromosome_2() == last_record.get_chromosome_2(), true, "The first and last records must have the same chromosome 2.");
            assert_eq!(first_record.get_strand_1() == last_record.get_strand_1(), true, "The first and last records must have the same strand 1.");
            assert_eq!(first_record.get_strand_2() == last_record.get_strand_2(), true, "The first and last records must have the same strand 2.");

            let graph_operation: GraphOperation = if *first_record.get_strand_1() == Strand::Forward {
                GraphOperation::new(
                    first_record.get_chromosome_1(),
                    first_record.get_position_1(),
                    first_record.get_strand_1().clone(),
                    GraphOperationType::Include,
                    last_record.get_chromosome_1(),
                    last_record.get_position_2(),
                    last_record.get_strand_1().clone(),
                    GraphOperationType::Include,
                    "".into(),
                    variant_type
                )
            } else {
                GraphOperation::new(
                    first_record.get_chromosome_1(),
                    first_record.get_position_2(),
                    first_record.get_strand_1().clone(),
                    GraphOperationType::Include,
                    last_record.get_chromosome_1(),
                    last_record.get_position_1(),
                    last_record.get_strand_1().clone(),
                    GraphOperationType::Include,
                    "".into(),
                    variant_type
                )
            };

            if included_read_bases.contains(&(first_record.get_start(), last_record.get_end())) == false {
                variant_records.push(
                    VariantRecord::new(
                        self.read_id,
                        first_record.get_start(),
                        last_record.get_end(),
                        graph_operation
                    )
                );
            }
        }

        variant_records
    }

    /// Identifies event variant records.
    ///
    /// This function identifies the following variant types:
    /// - Breakpoint
    /// - Circular RNA
    /// - Deletion
    /// - Exon truncation
    /// - Fusion gene
    /// - Translocation
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality.
    /// * `min_base_quality`: Minimum base quality.
    ///
    /// # Returns
    /// * Vector of VariantRecord objects.
    fn identify_event_variant_records(
        &self,
        records: &Vec<AlignmentStructureRecord>,
        min_mapping_quality: u16,
        min_base_quality: u8
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        for record in records.iter() {
            if *record.get_record_type() == AlignmentStructureRecordType::Base {
                continue;
            }

            // Check the mapping quality scores
            if record.get_mapping_quality_1() < min_mapping_quality || record.get_mapping_quality_2() < min_mapping_quality {
                continue;
            }

            match record.get_kind() {
                AlignmentStructureKind::Event(AlignmentStructureEventKind::Breakpoint) => {
                    // Get sequence
                    let sequence: &Box<str> = record.get_sequence();
                    let base_quality_scores: &Vec<u8> = record.get_base_quality_scores();
                    let mut filtered_sequence: String = String::new();
                    for (base, &quality) in sequence.chars().zip(base_quality_scores) {
                        if quality >= min_base_quality {
                            filtered_sequence.push(base);
                        }
                    }

                    // Get variant type
                    let mut variant_type: VariantType = VariantType::Breakpoint;
                    if record.get_context().is_some() {
                        match record.get_context().as_ref().unwrap() {
                            AlignmentStructureContext::Event(AlignmentStructureEventContext::BackSplicing) => {
                                variant_type = VariantType::CircularRNA;
                            },
                            AlignmentStructureContext::Event(AlignmentStructureEventContext::FusionGene) => {
                                variant_type = VariantType::FusionGene;
                            },
                            AlignmentStructureContext::Event(AlignmentStructureEventContext::NonCanonicalSplicing) => {
                                if record.get_chromosome_1() == record.get_chromosome_2() {
                                    variant_type = VariantType::Breakpoint;
                                } else {
                                    variant_type = VariantType::Translocation;
                                }
                            },
                            _ => {
                                // Do nothing
                            }
                        }
                    } else {
                        if record.get_chromosome_1() == record.get_chromosome_2() {
                            variant_type = VariantType::Breakpoint;
                        } else {
                            variant_type = VariantType::Translocation;
                        }
                    }

                    let graph_operation: GraphOperation = GraphOperation::new(
                        record.get_chromosome_1(),
                        record.get_position_1(),
                        record.get_strand_1().clone(),
                        record.get_operation_1().clone(),
                        record.get_chromosome_2(),
                        record.get_position_2(),
                        record.get_strand_2().clone(),
                        record.get_operation_2().clone(),
                        filtered_sequence.into(),
                        variant_type
                    );

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            record.get_start(),
                            record.get_end(),
                            graph_operation
                        )
                    );

                    if record.get_skipped().is_some() {
                        for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                            let reference_chromosome_id: u16 = reference_bases.first().unwrap().reference_chromosome_id;
                            let reference_position_1: u32 = reference_bases.first().unwrap().reference_position;
                            let reference_position_2: u32 = reference_bases.last().unwrap().reference_position;
                            let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                            let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                                reference_bases
                                    .iter()
                                    .map(|x| x.reference_nucleotide.as_str())
                                    .collect()
                            } else {
                                reference_bases
                                    .iter()
                                    .rev()
                                    .map(|x| x.reference_nucleotide.as_str())
                                    .collect()
                            };

                            let graph_operation: GraphOperation = GraphOperation::new(
                                reference_chromosome_id,
                                reference_position_1,
                                reference_strand.clone(),
                                GraphOperationType::Skip,
                                reference_chromosome_id,
                                reference_position_2,
                                reference_strand.clone(),
                                GraphOperationType::Skip,
                                sequence.into(),
                                VariantType::ExonTruncation
                            );

                            variant_records.push(
                                VariantRecord::new(
                                    self.read_id,
                                    record.get_start(),
                                    record.get_end(),
                                    graph_operation
                                )
                            );
                        }
                    }
                },
                AlignmentStructureKind::Event(AlignmentStructureEventKind::Boundary) => {
                    for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                        let reference_chromosome_id: u16 = reference_bases.first().unwrap().reference_chromosome_id;
                        let reference_position_1: u32 = reference_bases.first().unwrap().reference_position;
                        let reference_position_2: u32 = reference_bases.last().unwrap().reference_position;
                        let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                        let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                            reference_bases
                                .iter()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        } else {
                            reference_bases
                                .iter()
                                .rev()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        };

                        let graph_operation: GraphOperation = GraphOperation::new(
                            reference_chromosome_id,
                            reference_position_1,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            reference_chromosome_id,
                            reference_position_2,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            sequence.into(),
                            VariantType::ExonTruncation
                        );

                        variant_records.push(
                            VariantRecord::new(
                                self.read_id,
                                record.get_start(),
                                record.get_end(),
                                graph_operation
                            )
                        );
                    }
                },
                AlignmentStructureKind::Event(AlignmentStructureEventKind::Deletion) => {
                    let graph_operation: GraphOperation = GraphOperation::new(
                        record.get_chromosome_1(),
                        record.get_position_1(),
                        record.get_strand_1().clone(),
                        record.get_operation_1().clone(),
                        record.get_chromosome_2(),
                        record.get_position_2(),
                        record.get_strand_2().clone(),
                        record.get_operation_2().clone(),
                        "".into(),
                        VariantType::Deletion
                    );

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            record.get_start(),
                            record.get_end(),
                            graph_operation
                        )
                    );
                },
                AlignmentStructureKind::Event(AlignmentStructureEventKind::Splicing) => {
                    if record.get_context().is_some() {
                        if record.get_context().as_ref().unwrap() == &AlignmentStructureContext::Event(AlignmentStructureEventContext::FusionGene) {
                            // Get sequence
                            let sequence: &Box<str> = record.get_sequence();
                            let base_quality_scores: &Vec<u8> = record.get_base_quality_scores();
                            let mut filtered_sequence: String = String::new();
                            for (base, &quality) in sequence.chars().zip(base_quality_scores) {
                                if quality >= min_base_quality {
                                    filtered_sequence.push(base);
                                }
                            }

                            let graph_operation: GraphOperation = GraphOperation::new(
                                record.get_chromosome_1(),
                                record.get_position_1(),
                                record.get_strand_1().clone(),
                                record.get_operation_1().clone(),
                                record.get_chromosome_2(),
                                record.get_position_2(),
                                record.get_strand_2().clone(),
                                record.get_operation_2().clone(),
                                filtered_sequence.into(),
                                VariantType::FusionGene
                            );

                            variant_records.push(
                                VariantRecord::new(
                                    self.read_id,
                                    record.get_start(),
                                    record.get_end(),
                                    graph_operation
                                )
                            );
                        }
                    }

                    if record.get_skipped().is_some() {
                        for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                            let reference_chromosome_id: u16 = reference_bases.first().unwrap().reference_chromosome_id;
                            let reference_position_1: u32 = reference_bases.first().unwrap().reference_position;
                            let reference_position_2: u32 = reference_bases.last().unwrap().reference_position;
                            let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                            let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                                reference_bases
                                    .iter()
                                    .map(|x| x.reference_nucleotide.as_str())
                                    .collect()
                            } else {
                                reference_bases
                                    .iter()
                                    .rev()
                                    .map(|x| x.reference_nucleotide.as_str())
                                    .collect()
                            };

                            let graph_operation: GraphOperation = GraphOperation::new(
                                reference_chromosome_id,
                                reference_position_1,
                                reference_strand.clone(),
                                GraphOperationType::Skip,
                                reference_chromosome_id,
                                reference_position_2,
                                reference_strand.clone(),
                                GraphOperationType::Skip,
                                sequence.into(),
                                VariantType::ExonTruncation
                            );

                            variant_records.push(
                                VariantRecord::new(
                                    self.read_id,
                                    record.get_start(),
                                    record.get_end(),
                                    graph_operation
                                )
                            );
                        }
                    }
                },
                _ => {
                    continue;
                }
            };
        }

        variant_records
    }
}

impl Clone for AlignmentStructure {
    fn clone(&self) -> Self {
        AlignmentStructure {
            read_id: self.read_id,
            bases: self.bases.clone(),
            events: self.events.clone(),
            events_index: self.events_index.clone()
        }
    }
}
