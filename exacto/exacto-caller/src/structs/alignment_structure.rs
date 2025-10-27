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
    events: HashMap<(usize, usize), AlignmentStructureEvent>,

    /// Events index. The key is a read position and the value is another read position
    /// that together appear as a key in `events`.
    events_index: HashMap<usize, usize>
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
        // Index the event
        self.events_index
            .insert(event.get_prev_read_position(), event.get_next_read_position());
        self.events_index
            .insert(event.get_next_read_position(), event.get_prev_read_position());

        // Add the event
        self.events.insert(
            (event.get_prev_read_position(), event.get_next_read_position()),
            event
        );
    }

    pub fn contextualize(
        &mut self,
        reference_transcript_sequences: &Vec<&ReferenceTranscriptSequence>,
        gene_annotator: &(impl GeneAnnotator + Sync),
        chromosome_names_map: &BiMap<Box<str>, u16>
    ) {
        // Step 1. Make sure the reference transcript sequences have unique gene IDs
        let mut gene_ids: HashSet<Box<str>> = HashSet::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            if gene_ids.contains(reference_transcript_sequence.get_gene_id()) {
                panic!("Only 1 reference transcript per gene ID is allowed.");
            }
            gene_ids.insert(reference_transcript_sequence.get_gene_id().into());
        }

        // Step 2. Make sure each reference transcript sequence has 1 base that
        // overlaps with one of the self.bases
        let mut reference_transcript_overlap_map: HashMap<&str, bool> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            reference_transcript_overlap_map.insert(reference_transcript_sequence.get_transcript_id().into(), false);
        }
        if reference_transcript_sequences.is_empty() == false {
            for base in self.get_bases() {
                for reference_transcript_sequence in reference_transcript_sequences.iter() {
                    for base_reference in reference_transcript_sequence.get_bases() {
                        if base.get_reference_chromosome_id().unwrap() == base_reference.reference_chromosome_id &&
                            base.get_reference_position().unwrap() == base_reference.reference_position &&
                            base.get_reference_strand().as_ref().unwrap().clone() == base_reference.reference_strand {
                            reference_transcript_overlap_map.insert(reference_transcript_sequence.get_transcript_id(), true);
                        }
                    }
                }
            }

            // Make sure every reference transcript sequence overlaps with one of the self.bases
            for reference_transcript_sequence in reference_transcript_sequences.iter() {
                assert_eq!(
                    *reference_transcript_overlap_map.get(reference_transcript_sequence.get_transcript_id()).unwrap(), true,
                    "None of the ReferenceTranscriptSequence bases for {} overlaps with any of the AlignmentStructure bases.",
                    reference_transcript_sequence.get_transcript_id()
                );
            }
        }

        // Step 3. Index the positions of the reference transcript sequences
        let mut reference_transcripts_positions_map: HashMap<(u16, usize, Strand), &ReferenceTranscriptSequence> = HashMap::new();
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

        // Step 4. Identify context of each AlignmentStructureBase
        for i in 0..self.get_bases_length() {
            let base: &mut AlignmentStructureBase = self.get_mut_base(i);
            if *base.get_kind() != AlignmentStructureBaseKind::Unaligned {
                let key: (u16, usize, Strand) = (
                    base.get_reference_chromosome_id().unwrap(),
                    base.get_reference_position().unwrap(),
                    base.get_reference_strand().as_ref().unwrap().clone()
                );
                if reference_transcripts_positions_map.contains_key(&key) {
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
        let mut reference_transcripts_introns_map: HashMap<u16, HashSet<(usize, usize)>> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            let introns: Vec<(u16, usize, usize)> = reference_transcript_sequence.get_introns();
            for (chromosome_id, start, end) in introns.iter() {
                reference_transcripts_introns_map
                    .entry(*chromosome_id)
                    .or_insert_with(HashSet::new)
                    .insert((*start, *end));
            }
        }

        // Step 6. Identify context of each AlignmentStructureEvent
        let event_keys: Vec<(usize, usize)> = self.get_events().keys().cloned().collect();
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
                    let mut reference_start: usize = base_1.get_reference_position().unwrap();
                    let mut reference_end: usize = base_2.get_reference_position().unwrap();
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
                            self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::CanonicalSplicing);
                        } else {
                            if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                                self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                            } else {
                                self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                            }
                        }
                    } else {
                        self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                    }
                },
                AlignmentStructureEventKind::Breakpoint => {
                    if base_1.get_reference_chromosome_id().unwrap() == base_2.get_reference_chromosome_id().unwrap() {
                        let left_bases_set: HashSet<(u16, usize)> = (0..=read_position_1 as usize)
                            .map(|i| {
                                let base = self.get_base(i as usize);
                                (base.get_reference_chromosome_id().unwrap(), base.get_reference_position().unwrap())
                            })
                            .collect();
                        let right_bases_set: HashSet<(u16, usize)> = (read_position_2..self.get_bases_length())
                            .map(|i| {
                                let base = self.get_base(i);
                                (base.get_reference_chromosome_id().unwrap(), base.get_reference_position().unwrap())
                            })
                            .collect();
                        if left_bases_set.is_disjoint(&right_bases_set) == false {
                            self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::BackSplicing);
                        } else {
                            if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                                self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                            } else {
                                self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::NonCanonicalSplicing);
                            }
                        }
                    } else {
                        if gene_ids_disjoint_count_1 > 0 && gene_ids_disjoint_count_2 > 0 {
                            self.get_mut_event(read_position_1, read_position_2).set_context(AlignmentStructureEventContext::FusionGene);
                        }
                    }
                },
                _ => {
                    // Do nothing
                }
            }
        }

        // Step 7. Identify skipped reference transcript bases
        let mut reference_transcripts_bases: HashMap<(u16, usize, &Strand), (&str, &ReferenceBase)> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            for base in reference_transcript_sequence.get_bases() {
                reference_transcripts_bases.insert(
                    (base.reference_chromosome_id, base.reference_position, &base.reference_strand),
                    (reference_transcript_sequence.get_transcript_id(), base)
                );
            }
        }
        let alignment_structure_bases: HashSet<(u16, usize, &Strand)> = self
            .get_bases()
            .iter()
            .map(|base|
                (base.get_reference_chromosome_id().unwrap(),
                 base.get_reference_position().unwrap(),
                 base.get_reference_strand().as_ref().unwrap())
            )
            .collect();
        let mut reference_transcript_bases_skipped: HashSet<(u16, usize, &Strand)> = reference_transcripts_bases
            .keys()
            .cloned()
            .filter(|pos| !alignment_structure_bases.contains(pos))
            .sorted_by_key(|&(_, pos, _)| pos)
            .collect();

        // Step 8. Get the alignment structure base reference positions and sort them
        let mut base_reference_positions_map: HashMap<Box<str>, Vec<(usize, usize)>> = HashMap::new();
        for ((read_position_1, read_position_2), event) in self.get_events() {
            if event.get_kind() == &AlignmentStructureEventKind::Splicing {
                if *event.get_context().as_ref().unwrap() == AlignmentStructureEventContext::CanonicalSplicing {
                    continue;
                }
            }
            let base_1: &AlignmentStructureBase = self.get_base(*read_position_1);
            let base_2: &AlignmentStructureBase = self.get_base(*read_position_2);
            if base_1.get_reference_transcript_id().is_some() {
                let reference_transcript_id: Box<str> = base_1.get_reference_transcript_id().as_ref().unwrap().clone();
                base_reference_positions_map
                    .entry(reference_transcript_id)
                    .or_insert(Vec::new())
                    .push((*read_position_1, base_1.get_reference_position().unwrap()));
            }
            if base_2.get_reference_transcript_id().is_some() {
                let reference_transcript_id: Box<str> = base_2.get_reference_transcript_id().as_ref().unwrap().clone();
                base_reference_positions_map
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
                if base.is_embedded_insertion() == false && base.is_soft_clipped() == false {
                    base_reference_positions_map
                        .entry(reference_transcript_id.clone())
                        .or_insert(Vec::new())
                        .push(
                            (base.get_read_position(), base.get_reference_position().unwrap())
                        );
                    break;
                }
            }
            for base in reference_transcript_bases.iter().rev() {
                if base.is_embedded_insertion() == false && base.is_soft_clipped() == false {
                    base_reference_positions_map
                        .entry(reference_transcript_id.clone())
                        .or_insert(Vec::new())
                        .push(
                            (base.get_read_position(), base.get_reference_position().unwrap())
                        );
                    break;
                }
            }
        }
        for vec in base_reference_positions_map.values_mut() {
            vec.sort_by_key(|&(_, reference_position)| reference_position);
        }

        // Step 9. Identify the closet event for each skipped reference transcript base
        // Key: (chromosome ID, reference position, reference strand)
        // Value: (read position 1, read position 2)
        let mut events_map: HashMap<(u16, usize, &Strand), (usize, usize)> = HashMap::new();
        for (reference_chromosome_id, reference_position, reference_strand) in reference_transcript_bases_skipped.iter() {
            let reference_transcript_id: Box<str> = reference_transcripts_bases
                .get(&(*reference_chromosome_id, *reference_position, reference_strand))
                .unwrap()
                .0
                .into();

            // Identify the closest base
            let vec: &Vec<(usize, usize)> = base_reference_positions_map.get(&reference_transcript_id).unwrap();
            let closest_read_position: usize = match vec.binary_search_by_key(reference_position, |&(_, base_reference_position)| base_reference_position) {
                Ok(idx) => {
                    vec[idx].0
                },
                Err(idx) => {
                    if idx == 0 {
                        vec[0].0
                    } else if idx == vec.len() {
                        vec[vec.len() - 1].0
                    } else {
                        let read_position_1: usize = vec[idx-1].0;
                        let read_position_2: usize = vec[idx].0;
                        let reference_position_1: usize = vec[idx-1].1;
                        let reference_position_2: usize = vec[idx].1;
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
            } else {
                assert!(self.events_index.contains_key(&closest_read_position), "Event does not exist for base position {}. Skipped reference position: {}:{}", closest_read_position, reference_chromosome_id, reference_position);
                let closest_read_position_2: usize = *self.events_index.get(&closest_read_position).unwrap();
                assert!(self.has_event_between(closest_read_position, closest_read_position_2), "Event does not exist between bases {} and {}", closest_read_position, closest_read_position_2);
                events_map.insert(
                    (*reference_chromosome_id, *reference_position, reference_strand), (closest_read_position, closest_read_position_2)
                );
            }
        }

        assert_eq!(events_map.keys().len(), reference_transcript_bases_skipped.len(), "All skipped reference transcripts' bases should be mapped to an event");

        // Step 10. Insert skipped reference transcript bases into the alignment structure event
        for ((chromosome_id, reference_position, reference_strand), (read_position_1, read_position_2)) in events_map.iter() {
            let reference_base: &ReferenceBase = reference_transcripts_bases.get(&(*chromosome_id, *reference_position, reference_strand)).unwrap().1;
            assert!(self.has_event_between(*read_position_1, *read_position_2) == true);
            self.get_mut_event(*read_position_1, *read_position_2).add_skipped_reference_base(reference_base.clone());
        }
    }

    pub fn get_base(&self, read_position: usize) -> &AlignmentStructureBase {
        self.bases.get(read_position as usize).unwrap()
    }

    pub fn get_bases(&self) -> &Vec<AlignmentStructureBase> {
        &self.bases
    }

    pub fn get_bases_length(&self) -> usize {
        self.bases.len() as usize
    }
    
    pub fn get_event(&self,
        read_position_1: usize,
        read_position_2: usize
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

    pub fn get_events(&self) -> &HashMap<(usize, usize), AlignmentStructureEvent> {
        &self.events
    }

    pub fn get_event_at_read_position(&self, read_position: usize) -> &AlignmentStructureEvent {
        let read_position_2 = self.events_index.get(&read_position).unwrap();
        self.get_event(read_position, *read_position_2)
    }
    
    pub fn get_events_of_kind(&self, kind: AlignmentStructureEventKind) -> Vec<(usize, usize, &AlignmentStructureEvent)> {
        let mut events: Vec<(usize, usize, &AlignmentStructureEvent)> = Vec::new();
        for ((read_position_1, read_position_2), event) in self.get_events() {
            if event.get_kind() == &kind {
                events.push((*read_position_1, *read_position_2, event));
            }
        }
        events
    }

    pub fn get_mut_base(&mut self, read_position: usize) -> &mut AlignmentStructureBase {
        self.bases.get_mut(read_position as usize).unwrap()
    }

    pub fn get_mut_event(
        &mut self,
        read_position_1: usize,
        read_position_2: usize
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

    pub fn has_event(&self, read_position: usize) -> bool {
        if self.events_index.contains_key(&read_position) {
            true
        } else {
            false
        }
    }

    pub fn has_event_between(
        &self,
        read_position_1: usize,
        read_position_2: usize
    ) -> bool {
        if read_position_1 < read_position_2 {
            self.events.contains_key(&(read_position_1, read_position_2))
        } else {
            self.events.contains_key(&(read_position_2, read_position_1))
        }
    }

    pub fn identify_exons(&self) -> Vec<TranscriptModelExon> {
        // Step 1. Cluster adjacent bases by reference position
        let mut uf: UnionFind = UnionFind::new();
        for i in 0..self.get_bases_length() {
            uf.union(i as usize, i as usize);
            if i > 0 {
                let prev_base: &AlignmentStructureBase = self.get_base(i - 1);
                let curr_base: &AlignmentStructureBase = self.get_base(i);
                if self.has_event_between(i-1, i) {
                    let event: &AlignmentStructureEvent = self.get_event(i-1, i);
                    if event.get_kind() == &AlignmentStructureEventKind::Deletion {
                        if prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                            prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() {
                            uf.union(i as usize - 1, i as usize);
                        }
                    }
                } else {
                    if prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                        prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() &&
                        prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1 {
                        uf.union(i as usize - 1, i as usize);
                    }
                }
            }
        }

        // Step 2. Identify exonic boundaries
        let mut exons: Vec<TranscriptModelExon> = Vec::new();
        let mut exon_number: u16 = 1;
        let mut clusters: Vec<HashSet<usize>> = uf.get_clusters();
        for cluster in clusters.iter() {
            // Sort the read positions
            let mut read_positions: Vec<usize> = cluster.iter().map(|&pos| pos as usize).collect();
            read_positions.sort();

            let mut bases: Vec<&AlignmentStructureBase> = Vec::new();
            for read_position in read_positions.iter() {
                bases.push(self.get_base(*read_position));
            }

            let reference_chromosome_id: u16 = bases.first().unwrap().get_reference_chromosome_id().unwrap();
            let reference_strand: Strand = bases.first().unwrap().get_reference_strand().as_ref().unwrap().clone();
            let reference_start: usize = bases.iter().map(|base| base.get_reference_position().unwrap()).min().unwrap();
            let reference_end: usize = bases.iter().map(|base| base.get_reference_position().unwrap()).max().unwrap();
            let read_start_position: usize = bases.first().unwrap().get_read_position();
            let read_end_position: usize = bases.last().unwrap().get_read_position();

            let exon: TranscriptModelExon = TranscriptModelExon::new(
                reference_chromosome_id,
                reference_start,
                reference_end,
                reference_strand,
                exon_number,
                read_start_position,
                read_end_position
            );

            exons.push(exon);
            exon_number += 1;
        }

        exons
    }

    pub fn identify_introns(
        &self,
        chromosome_names_map: &BiMap<Box<str>,u16>,
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
                        reference_start_position as usize + 1,
                        reference_start_position as usize + 2,
                        reference_genome_fasta_file
                    );

                    let sequence_2: Box<str> = get_fasta_sequence(
                        &*reference_chromosome_name,
                        reference_end_position as usize - 2,
                        reference_end_position as usize - 1,
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
        // Step 1. Cluster bases
        let num_bases: usize = self.get_bases_length();
        let mut uf_bases: UnionFind = UnionFind::new();
        for i in 0..num_bases {
            let curr_base: &AlignmentStructureBase = self.get_base(i);
            if !curr_base.is_embedded_insertion() {
                uf_bases.union(i as usize, i as usize);
            }
        }
        for i in 0..num_bases {
            if i > 0 {
                let prev_base: &AlignmentStructureBase = self.get_base(i - 1);
                let curr_base: &AlignmentStructureBase = self.get_base(i);
                if !prev_base.is_embedded_insertion() && !curr_base.is_embedded_insertion() {
                    if prev_base.get_context().is_some() && curr_base.get_context().is_some() {
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            *prev_base.get_context().as_ref().unwrap() == *curr_base.get_context().as_ref().unwrap() &&
                            prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                            prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1 &&
                            prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() {
                            uf_bases.union(i as usize - 1, i as usize);
                        }
                    } else {
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            prev_base.get_reference_chromosome_id().unwrap() == curr_base.get_reference_chromosome_id().unwrap() &&
                            prev_base.get_reference_position().unwrap().abs_diff(curr_base.get_reference_position().unwrap()) <= 1 &&
                            prev_base.get_reference_strand().as_ref().unwrap() == curr_base.get_reference_strand().as_ref().unwrap() {
                            uf_bases.union(i as usize - 1, i as usize);
                        }
                    }
                }
            }
        }

        // Step 2. Record bases
        let mut records: Vec<AlignmentStructureRecord> = Vec::new();
        for cluster in uf_bases.get_clusters() {
            let mut read_positions: Vec<usize> = cluster.into_iter().collect();
            read_positions.sort();
            let bases: Vec<&AlignmentStructureBase> = read_positions.iter().map(|&i| self.get_base(i as usize)).collect();
            let first_base: &AlignmentStructureBase = bases.first().unwrap();
            let last_base: &AlignmentStructureBase = bases.last().unwrap();
            let mut sequence: String = String::new();
            let mut base_quality_scores: Vec<u8> = Vec::new();
            for base in bases.iter() {
                sequence.push_str(base.get_nucleotide().as_str());
                base_quality_scores.push(base.get_base_quality());
            }
            let (base_1, base_2) = if *first_base.get_kind() == AlignmentStructureBaseKind::Unaligned {
                (first_base, last_base)
            } else {
                if *first_base.get_reference_strand().as_ref().unwrap() == Strand::Forward {
                    (first_base, last_base)
                } else {
                    (last_base, first_base)
                }
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
                    if first_base.is_embedded_insertion() == false && last_base.is_embedded_insertion() == false {
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
                    }
                },
                AlignmentStructureBaseKind::Unaligned => {
                    let record: AlignmentStructureRecord = AlignmentStructureRecord::new(
                        first_base.get_read_position(),
                        last_base.get_read_position(),
                        sequence.as_str(),
                        base_quality_scores,
                        AlignmentStructureRecordType::Base,
                        AlignmentStructureKind::Base(first_base.get_kind().clone()),
                        None,
                        0,
                        0,
                        GraphOperationType::Noop,
                        Strand::Unknown,
                        0,
                        0,
                        0,
                        GraphOperationType::Noop,
                        Strand::Unknown,
                        0,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None,
                        None
                    );
                    records.push(record);
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
        min_mapping_quality: usize,
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

        // Step 3. Identify event variant records
        variant_records.extend(self.identify_event_variant_records(
            &records,
            min_mapping_quality,
            min_base_quality
        ));

        // Step 4. Sort the variant records
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
        min_mapping_quality: usize,
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

                    let decrement: usize = sequence.len().abs_diff(curr_record.get_sequence().len()) as usize;

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
                    // Exclude template-switching artifacts if the analyte type is RNA
                    if analyte_type == AnalyteType::RNA &&
                        (curr_record.get_start() == 0 || curr_record.get_end() == self.get_bases_length() - 1) &&
                        curr_record.get_sequence().len() <= 3 &&
                        (curr_record.get_sequence().to_uppercase().contains("AA") ||
                         curr_record.get_sequence().to_uppercase().contains("CC") ||
                         curr_record.get_sequence().to_uppercase().contains("GG") ||
                         curr_record.get_sequence().to_uppercase().contains("TT")) {
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

        if analyte_type == AnalyteType::DNA {
            return variant_records;
        }

        // Step 2. Identify variant records based on the AlignmentStructureBase contexts
        let mut uf: UnionFind = UnionFind::new();
        for i in 0..records.len() {
            let record: &AlignmentStructureRecord = records.get(i).unwrap();
            if *record.get_record_type() == AlignmentStructureRecordType::Base &&
                *record.get_context().as_ref().unwrap() != AlignmentStructureContext::Base(AlignmentStructureBaseContext::Exonic) {
                uf.union(i, i);
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
                            uf.union(i - 1, i);
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
                            uf.union(i - 1, i);
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
                        uf.union(i - 1, i);
                    }
                }
            }
        }
        for cluster in uf.get_clusters().iter() {
            // Sort the record positions
            let mut record_indices: Vec<usize> = cluster.iter().map(|&pos| pos as usize).collect();
            record_indices.sort();

            let first_record: &AlignmentStructureRecord = records.get(*record_indices.first().unwrap() as usize).unwrap();
            let last_record: &AlignmentStructureRecord = records.get(*record_indices.last().unwrap() as usize).unwrap();

            let prev_record: Option<&AlignmentStructureRecord> = if *record_indices.first().unwrap() > 0usize  {
                records.get(record_indices[0] as usize - 1)
            } else {
                None
            };

            let next_record: Option<&AlignmentStructureRecord> = if *record_indices.last().unwrap() < records.len() as usize - 1  {
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

            variant_records.push(
                VariantRecord::new(
                    self.read_id,
                    first_record.get_start(),
                    last_record.get_end(),
                    graph_operation
                )
            );
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
        min_mapping_quality: usize,
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
                    let sequence = record.get_sequence();
                    let base_quality_scores = record.get_base_quality_scores();
                    let mut filtered_sequence = String::new();
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
                            let reference_position_1: usize = reference_bases.first().unwrap().reference_position;
                            let reference_position_2: usize = reference_bases.last().unwrap().reference_position;
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
                        let reference_position_1: usize = reference_bases.first().unwrap().reference_position;
                        let reference_position_2: usize = reference_bases.last().unwrap().reference_position;
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
                            let sequence = record.get_sequence();
                            let base_quality_scores = record.get_base_quality_scores();
                            let mut filtered_sequence = String::new();
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
                            let reference_position_1: usize = reference_bases.first().unwrap().reference_position;
                            let reference_position_2: usize = reference_bases.last().unwrap().reference_position;
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
