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
        let mut overlaps: HashMap<&str, bool> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            overlaps.insert(reference_transcript_sequence.get_transcript_id().into(), false);
        }
        if reference_transcript_sequences.is_empty() == false {
            for base in self.get_bases() {
                for reference_transcript_sequence in reference_transcript_sequences.iter() {
                    for base_reference in reference_transcript_sequence.get_bases() {
                        if base.get_reference_chromosome_id().unwrap() == base_reference.reference_chromosome_id &&
                            base.get_reference_position().unwrap() == base_reference.reference_position &&
                            base.get_reference_strand().as_ref().unwrap().clone() == base_reference.reference_strand {
                            overlaps.insert(reference_transcript_sequence.get_transcript_id(), true);
                        }
                    }
                }
            }

            // Make sure every reference transcript sequence overlaps with one of the self.bases
            for reference_transcript_sequence in reference_transcript_sequences.iter() {
                assert_eq!(
                    *overlaps.get(reference_transcript_sequence.get_transcript_id()).unwrap(), true,
                    "None of the ReferenceTranscriptSequence bases for {} overlaps with any of the AlignmentStructure bases.",
                    reference_transcript_sequence.get_transcript_id()
                );
            }
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

        // Step 4. Identify context of each AlignmentStructureBase
        for i in 0..self.get_bases_length() {
            let base: &mut AlignmentStructureBase = self.get_mut_base(i);
            if *base.get_kind() != AlignmentStructureBaseKind::Unaligned {
                let key: (u16, u32, Strand) = (
                    base.get_reference_chromosome_id().unwrap(),
                    base.get_reference_position().unwrap(),
                    base.get_reference_strand().as_ref().unwrap().clone()
                );
                if reference_transcripts_positions_map.contains_key(&key) {
                    let reference_transcript_sequence: &ReferenceTranscriptSequence = reference_transcripts_positions_map.get(&key).unwrap();
                    base.set_context(AlignmentStructureBaseContext::Exonic);
                    base.set_reference_transcript_id(reference_transcript_sequence.get_transcript_id());
                } else {
                    for reference_transcript_sequence in reference_transcript_sequences.iter() {
                        if base.get_reference_position().unwrap() >= reference_transcript_sequence.get_transcript_start() &&
                            base.get_reference_position().unwrap() <= reference_transcript_sequence.get_transcript_end() {
                            base.set_context(AlignmentStructureBaseContext::Intronic);
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

        // Step 6. Identify context of each AlignmentStructureEvent
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
                        let left_bases_set: HashSet<(u16, u32)> = (0..=read_position_1 as usize)
                            .map(|i| {
                                let base = self.get_base(i as u32);
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
        let mut reference_transcripts_bases: HashMap<(u16, u32, &Strand), (&str, &ReferenceBase)> = HashMap::new();
        for reference_transcript_sequence in reference_transcript_sequences.iter() {
            for base in reference_transcript_sequence.get_bases() {
                reference_transcripts_bases.insert(
                    (base.reference_chromosome_id, base.reference_position, &base.reference_strand),
                    (reference_transcript_sequence.get_transcript_id(), base)
                );
            }
        }
        let alignment_structure_bases: HashSet<(u16, u32, &Strand)> = self
            .get_bases()
            .iter()
            .map(|base|
                (base.get_reference_chromosome_id().unwrap(),
                 base.get_reference_position().unwrap(),
                 base.get_reference_strand().as_ref().unwrap())
            )
            .collect();
        let mut reference_transcript_bases_skipped: HashSet<(u16, u32, &Strand)> = reference_transcripts_bases
            .keys()
            .cloned()
            .filter(|pos| !alignment_structure_bases.contains(pos))
            .sorted_by_key(|&(_, pos, _)| pos)
            .collect();

        // Step 8. Get the alignment structure base reference positions and sort them
        let mut base_reference_positions_map: HashMap<Box<str>, Vec<(u32, u32)>> = HashMap::new();
        for ((read_position_1, read_position_2), event) in self.get_events() {
            if event.get_kind() == &AlignmentStructureEventKind::Splicing {
                if *event.get_context().unwrap() == AlignmentStructureEventContext::CanonicalSplicing {
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
        let mut events_map: HashMap<(u16, u32, &Strand), (u32, u32)> = HashMap::new();
        for (reference_chromosome_id, reference_position, reference_strand) in reference_transcript_bases_skipped.iter() {
            let reference_transcript_id: Box<str> = reference_transcripts_bases
                .get(&(*reference_chromosome_id, *reference_position, reference_strand))
                .unwrap()
                .0
                .into();

            // Identify the closest base
            let vec: &Vec<(u32, u32)> = base_reference_positions_map.get(&reference_transcript_id).unwrap();
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
            } else {
                assert!(self.events_index.contains_key(&closest_read_position), "Event does not exist for base position {}. Skipped reference position: {}:{}", closest_read_position, reference_chromosome_id, reference_position);
                let closest_read_position_2: u32 = *self.events_index.get(&closest_read_position).unwrap();
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
            self.get_mut_event(*read_position_1, *read_position_2).add_reference_base(reference_base.clone());
        }
    }

    pub fn get_base(&self, read_position: u32) -> &AlignmentStructureBase {
        self.bases.get(read_position as usize).unwrap()
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

    pub fn get_mut_base(&mut self, read_position: u32) -> &mut AlignmentStructureBase {
        self.bases.get_mut(read_position as usize).unwrap()
    }

    pub fn get_mut_event(
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

    /// Identifies DNA variant records.
    ///
    /// The following DNA variant types are identified:
    ///     - Single-nucleotide variant
    ///     - Multi-nucleotide variant
    ///     - Insertion
    ///     - Deletion
    ///     - Breakpoint
    ///     - Translocation
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality (inclusive).
    /// * `min_base_quality`: Minimum base quality (inclusive).
    ///
    /// # Returns
    /// Vector of VariantRecord objects.
    pub fn identify_dna_variant_records(
        &self,
        min_mapping_quality: u32,
        min_base_quality: u8
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = self.identify_base_kind_variant_records(
            min_mapping_quality,
            min_base_quality
        );
        variant_records.extend(self.identify_event_kind_variant_records(min_mapping_quality));
        variant_records.sort_by(|a, b| {
            a.get_chromosome_1()
                .cmp(&b.get_chromosome_1())
                .then(a.get_position_1().cmp(&b.get_position_1()))
        });
        variant_records
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
            let mut read_positions: Vec<u32> = cluster.iter().map(|&pos| pos as u32).collect();
            read_positions.sort();

            let mut bases: Vec<&AlignmentStructureBase> = Vec::new();
            for read_position in read_positions.iter() {
                bases.push(self.get_base(*read_position));
            }

            let reference_chromosome_id: u16 = bases.first().unwrap().get_reference_chromosome_id().unwrap();
            let reference_strand: Strand = bases.first().unwrap().get_reference_strand().as_ref().unwrap().clone();
            let reference_start: u32 = bases.iter().map(|base| base.get_reference_position().unwrap()).min().unwrap();
            let reference_end: u32 = bases.iter().map(|base| base.get_reference_position().unwrap()).max().unwrap();
            let read_start_position: u32 = bases.first().unwrap().get_read_position();
            let read_end_position: u32 = bases.last().unwrap().get_read_position();

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

    /// Identifies RNA variant records.
    ///
    /// The following RNA variant types are identified:
    ///     - Single-nucleotide variant
    ///     - Multi-nucleotide variant
    ///     - Insertion
    ///     - Deletion
    ///     - Breakpoint
    ///     - Translocation
    ///     - Fusion gene
    ///     - Circular RNA
    ///     - Cryptic exon
    ///     - Exon truncation
    ///     - Intron retention
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality (inclusive).
    /// * `min_base_quality`: Minimum base quality (inclusive).
    /// * `gene_annotator`: Gene annotator.
    /// * `reference_transcript_sequence`: Reference to ReferenceTranscriptSequence object.
    /// * `chromosome_names_map`: Chromosome names bimap.
    ///
    /// # Returns
    /// Vector of VariantRecord objects.
    pub fn identify_rna_variant_records(
        &self,
        min_mapping_quality: u32,
        min_base_quality: u8,
        gene_annotator: &(impl GeneAnnotator + Sync),
    ) -> Vec<VariantRecord> {
        // Step 1. Identify variant records based on the AlignmentStructureBase kinds
        let mut variant_records: Vec<VariantRecord> = self.identify_base_kind_variant_records(
            min_mapping_quality,
            min_base_quality
        );

        // Step 2. Identify variant records based on the AlignmentStructureBase contexts
        variant_records.extend(self.identify_base_context_variant_records(
            min_mapping_quality,
            min_base_quality
        ));

        // Step 3. Identify variant records based on the AlignmentStructureEvent kinds
        variant_records.extend(self.identify_event_kind_variant_records(min_mapping_quality));

        // Step 4. Identify variant records based on the AlignmentStructureEvent contexts
        variant_records.extend(self.identify_event_context_variant_records(
            min_mapping_quality,
            gene_annotator
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

    // pub fn to_string(&self) -> String {
    //     let mut structure: String = "".to_string();
    //     for i in 0..self.get_bases_length() {
    //         let base: &AlignmentStructureBase = self.get_base(i);
    //     }
    // }
}

/// Helper functions
impl AlignmentStructure {

    /// Identifies variant records in self.bases based on the AlignmentStructureBaseKind.
    ///
    /// This function identifies the following variant types:
    /// - Single-nucleotide variant
    /// - Multi-nucleotide variant
    /// - Insertion.
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality.
    /// * `min_base_quality`: Minimum base quality.
    ///
    /// # Returns
    /// * Vector of VariantRecord objects.
    fn identify_base_kind_variant_records(
        &self,
        min_mapping_quality: u32,
        min_base_quality: u8
    ) -> Vec<VariantRecord> {
        // Step 1. Cluster the bases based on the base kind
        let mut uf: UnionFind = UnionFind::new();
        for base in self.get_bases().iter() {
            let id = base.get_read_position() as usize;
            match *base.get_kind() {
                AlignmentStructureBaseKind::Mismatch => {
                    uf.union(id, id);
                },
                AlignmentStructureBaseKind::Insertion => {
                    if base.is_embedded_insertion() == false {
                        uf.union(id, id);
                    }
                },
                _ => {
                    // Do nothing
                }
            }
        }
        for pair in self.get_bases().windows(2) {
            // Union adjacent bases that are contiguous and share the same context
            let prev_base: &AlignmentStructureBase = &pair[0];
            let curr_base: &AlignmentStructureBase = &pair[1];

            if prev_base.is_embedded_insertion() || curr_base.is_embedded_insertion() {
                continue;
            }

            if curr_base.get_mapping_quality().unwrap() < min_mapping_quality {
                continue;
            }

            if curr_base.get_base_quality() < min_base_quality {
                continue;
            }

            // Make sure the current base is not an embedded insertion
            if curr_base.is_embedded_insertion() {
                continue;
            }

            // Require all reference fields to be present on both sides
            if let (
                Some(&prev_chr),
                Some(&curr_chr),
                Some(prev_strand),
                Some(curr_strand),
                Some(prev_pos),
                Some(curr_pos),
            ) = (
                prev_base.get_reference_chromosome_id().as_ref(),
                curr_base.get_reference_chromosome_id().as_ref(),
                prev_base.get_reference_strand().as_ref(),
                curr_base.get_reference_strand().as_ref(),
                prev_base.get_reference_position().as_ref(),
                curr_base.get_reference_position().as_ref(),
            ) {
                if prev_chr == curr_chr && prev_strand == curr_strand {
                    // MNV
                    if prev_pos.abs_diff(*curr_pos) == 1 {
                        if *prev_base.get_kind() == *curr_base.get_kind() && *prev_base.get_kind() == AlignmentStructureBaseKind::Mismatch {
                            uf.union(
                                prev_base.get_read_position() as usize,
                                curr_base.get_read_position() as usize,
                            );
                        }
                    }

                    // Insertion
                    if *prev_pos == *curr_pos {
                        if *prev_base.get_kind() == *curr_base.get_kind() && *prev_base.get_kind() == AlignmentStructureBaseKind::Insertion {
                            uf.union(
                                prev_base.get_read_position() as usize,
                                curr_base.get_read_position() as usize,
                            );
                        }
                    }
                }
            }
        }

        // Step 2. Identify SNVs, MNVs, and insertions
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        let read_sequence_length: u32 = self.get_bases_length();
        for cluster in uf.get_clusters().iter() {
            // Sort the read positions
            let mut read_positions: Vec<u32> = cluster.iter().map(|&pos| pos as u32).collect();
            read_positions.sort();

            let bases: Vec<&AlignmentStructureBase> = read_positions.iter().map(|&pos| self.get_base(pos)).collect();
            if read_positions.len() == 1 {
                assert!(*bases[0].get_kind() == AlignmentStructureBaseKind::Mismatch);

                // Exclude SNVs at the first or last base of the read
                if bases[0].get_read_position() == 0 ||
                    bases[0].get_read_position() == read_sequence_length - 1 {
                    continue;
                }

                let graph_operation: GraphOperation = GraphOperation::new(
                    bases[0].get_reference_chromosome_id().unwrap(),
                    bases[0].get_reference_position().unwrap() - 1,
                    bases[0].get_reference_strand().as_ref().unwrap().clone(),
                    GraphOperationType::Downstream,
                    bases[0].get_reference_chromosome_id().unwrap(),
                    bases[0].get_reference_position().unwrap() + 1,
                    bases[0].get_reference_strand().as_ref().unwrap().clone(),
                    GraphOperationType::Upstream,
                    bases[0].get_nucleotide().as_str().into(),
                    VariantType::SingleNucleotideVariant
                );
                variant_records.push(
                    VariantRecord::new(
                        self.read_id,
                        bases[0].get_read_position(),
                        bases[0].get_read_position(),
                        graph_operation
                    )
                );
            } else {
                if *bases[0].get_kind() == AlignmentStructureBaseKind::Mismatch {
                    let first_base: &AlignmentStructureBase = bases.first().unwrap();
                    let last_base: &AlignmentStructureBase = bases.last().unwrap();
                    let sequence: String = bases
                        .iter()
                        .map(|b| b.get_nucleotide().as_str())
                        .collect::<Vec<_>>()
                        .join("");
                    let (chrom1, pos1, strand1, op1, chrom2, pos2, strand2, op2, seq_box) = if bases[0].get_reference_strand().as_ref().unwrap().clone() == Strand::Forward {
                        (
                            first_base.get_reference_chromosome_id().unwrap(),
                            first_base.get_reference_position().unwrap() - 1,
                            first_base.get_reference_strand().as_ref().unwrap().clone(),
                            GraphOperationType::Downstream,
                            last_base.get_reference_chromosome_id().unwrap(),
                            last_base.get_reference_position().unwrap() + 1,
                            last_base.get_reference_strand().as_ref().unwrap().clone(),
                            GraphOperationType::Upstream,
                            sequence.into_boxed_str()
                        )
                    } else {
                        (
                            last_base.get_reference_chromosome_id().unwrap(),
                            last_base.get_reference_position().unwrap() + 1,
                            last_base.get_reference_strand().as_ref().unwrap().clone(),
                            GraphOperationType::Downstream,
                            first_base.get_reference_chromosome_id().unwrap(),
                            first_base.get_reference_position().unwrap() - 1,
                            first_base.get_reference_strand().as_ref().unwrap().clone(),
                            GraphOperationType::Upstream,
                            reverse_complement(&sequence)
                        )
                    };
                    let graph_operation: GraphOperation = GraphOperation::new(
                        chrom1, pos1, strand1, op1,
                        chrom2, pos2, strand2, op2,
                        seq_box,
                        VariantType::MultiNucleotideVariant
                    );
                    variant_records.push(VariantRecord::new(
                        self.read_id,
                        first_base.get_read_position(),
                        last_base.get_read_position(),
                        graph_operation,
                    ));
                }
                if *bases[0].get_kind() == AlignmentStructureBaseKind::Insertion {
                    let mut sequence: String = bases
                        .iter()
                        .map(|b| b.get_nucleotide().as_str())
                        .collect::<Vec<_>>()
                        .join("");
                    let first_base: &AlignmentStructureBase = bases.first().unwrap();
                    let last_base: &AlignmentStructureBase = bases.last().unwrap();
                    let is_softclipped = bases.iter().any(|b| b.is_soft_clipped());
                    if is_softclipped {
                        assert!(
                            first_base.get_read_position() == 0 ||
                            last_base.get_read_position() == read_sequence_length - 1,
                            "Softclipped bases are expected to be either at \
                            the start or end of the read sequence."
                        );

                        // Skip likely template-switching artifacts
                        if sequence.len() == 1 {
                            continue;
                        }
                        if sequence.len() <= 3 {
                            let up: String = sequence.to_uppercase();
                            if up.contains("AA") ||
                                up.contains("CC") ||
                                up.contains("GG") ||
                                up.contains("TT") {
                                continue;
                            }
                        }
                        let chromosome_id: u16 = first_base.get_reference_chromosome_id().unwrap();
                        let strand: Strand = first_base.get_reference_strand().as_ref().unwrap().clone();
                        let at_start: bool = first_base.get_read_position() == 0;
                        let (reference_position_1, reference_position_2) = match (at_start, strand.clone()) {
                            (true, Strand::Forward)  => (first_base.get_reference_position().unwrap() - 1, first_base.get_reference_position().unwrap()),
                            (true, Strand::Reverse)  => (first_base.get_reference_position().unwrap(),     first_base.get_reference_position().unwrap() + 1),
                            (false, Strand::Forward) => (first_base.get_reference_position().unwrap(),     first_base.get_reference_position().unwrap() + 1),
                            (false, Strand::Reverse) => (first_base.get_reference_position().unwrap() - 1, first_base.get_reference_position().unwrap()),
                            _ => panic!("Invalid strand")
                        };
                        let graph_operation = GraphOperation::new(
                            chromosome_id, reference_position_1, strand.clone(), GraphOperationType::Downstream,
                            chromosome_id, reference_position_2, strand.clone(), GraphOperationType::Upstream,
                            sequence.into_boxed_str(),
                            VariantType::Insertion
                        );
                        variant_records.push(VariantRecord::new(
                            self.read_id,
                            first_base.get_read_position(),
                            last_base.get_read_position(),
                            graph_operation,
                        ));
                    } else {
                        let chromosome_id: u16 = first_base.get_reference_chromosome_id().unwrap();
                        let reference_position: u32 = first_base.get_reference_position().unwrap();
                        let strand: Strand = first_base.get_reference_strand().as_ref().unwrap().clone();
                        let graph_operation = GraphOperation::new(
                            chromosome_id, reference_position, strand.clone(), GraphOperationType::Downstream,
                            chromosome_id, reference_position + 1, strand.clone(), GraphOperationType::Upstream,
                            sequence.into_boxed_str(),
                            VariantType::Insertion,
                        );
                        variant_records.push(VariantRecord::new(
                            self.read_id,
                            first_base.get_read_position(),
                            last_base.get_read_position(),
                            graph_operation
                        ));
                    }
                }
            }
        }

        variant_records
    }

    /// Identifies variant records in self.bases based on the AlignmentStructureBaseContext.
    ///
    /// This function identifies the following variant types:
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
    fn identify_base_context_variant_records(
        &self,
        min_mapping_quality: u32,
        min_base_quality: u8
    ) -> Vec<VariantRecord> {
        // Step 1. Cluster the bases based on the base context
        let mut uf: UnionFind = UnionFind::new();
        for base in self.get_bases().iter() {
            if base.get_context().is_some() {
                if *base.get_context().as_ref().unwrap() != AlignmentStructureBaseContext::Exonic {
                    let id = base.get_read_position() as usize;
                    uf.union(id, id);
                }
            }
        }
        for pair in self.get_bases().windows(2) {
            // Union adjacent bases that are contiguous and share the same context
            let prev: &AlignmentStructureBase = &pair[0];
            let curr: &AlignmentStructureBase = &pair[1];

            if curr.get_mapping_quality().unwrap() < min_mapping_quality {
                continue;
            }

            if curr.get_base_quality() < min_base_quality {
               continue;
            }

            // Require all reference fields to be present on both sides
            if let (
                Some(&prev_chr),
                Some(&curr_chr),
                Some(prev_strand),
                Some(curr_strand),
                Some(prev_pos),
                Some(curr_pos),
            ) = (
                prev.get_reference_chromosome_id().as_ref(),
                curr.get_reference_chromosome_id().as_ref(),
                prev.get_reference_strand().as_ref(),
                curr.get_reference_strand().as_ref(),
                prev.get_reference_position().as_ref(),
                curr.get_reference_position().as_ref(),
            ) {
                if prev_chr == curr_chr &&
                    prev_strand == curr_strand &&
                    prev_pos.abs_diff(*curr_pos) == 1 &&
                    *prev.get_context().as_ref().unwrap() == *curr.get_context().as_ref().unwrap() &&
                    prev.get_context().as_ref().unwrap().clone() != AlignmentStructureBaseContext::Exonic {
                    uf.union(
                        prev.get_read_position() as usize,
                        curr.get_read_position() as usize
                    );
                }
            }
        }

        // Step 2. Identify cryptic exons and intron retentions
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        for cluster in uf.get_clusters().iter() {
            // Sort the read positions
            let mut read_positions: Vec<u32> = cluster.iter().map(|&pos| pos as u32).collect();
            read_positions.sort();

            // Get the first and last bases in the cluster
            let first_base: &AlignmentStructureBase = self.get_base(read_positions[0]);
            let last_base: &AlignmentStructureBase = self.get_base(read_positions[read_positions.len() - 1]);

            // Get the previous and next bases
            let prev_base: Option<&AlignmentStructureBase> = if first_base.get_read_position() > 0 {
                Some(self.get_base(first_base.get_read_position() - 1))
            } else {
                None
            };

            let next_base: Option<&AlignmentStructureBase> = if last_base.get_read_position() < self.get_bases_length() - 1 {
                Some(self.get_base(last_base.get_read_position() + 1))
            } else {
                None
            };

            // Determine if the cluster is a cryptic exon, intron retention, or UTR extension
            // If there is an exonic base immediately adjacent in reference coordinates
            // (same chromosome, same strand, and abs_diff(ref_pos) <= 1), classify this
            // cluster of read positions as an intron retention.
            // If there is an intergenic base immediately adjacent in reference coordinates
            // (same chromosome, same strand, and abs_diff(ref_pos) <= 1), classify this
            // cluster of read positions as an UTR extension.
            // Otherwise, classify as a cryptic exon.
            let mut variant_type: VariantType = VariantType::CrypticExon;
            if prev_base.is_some() {
                let prev_base_: &AlignmentStructureBase = prev_base.unwrap();
                if prev_base_.get_reference_chromosome_id().unwrap() == first_base.get_reference_chromosome_id().unwrap() &&
                    *prev_base_.get_reference_strand().as_ref().unwrap() == *first_base.get_reference_strand().as_ref().unwrap() &&
                    prev_base_.get_reference_position().unwrap().abs_diff(first_base.get_reference_position().unwrap()) <= 1 &&
                    prev_base_.get_context().as_ref().unwrap().clone() == AlignmentStructureBaseContext::Exonic {
                    if first_base.get_context().as_ref().unwrap() == &AlignmentStructureBaseContext::Intronic {
                        variant_type = VariantType::IntronRetention;
                    }
                    if first_base.get_context().as_ref().unwrap() == &AlignmentStructureBaseContext::Intergenic {
                        variant_type = VariantType::UTRExtension;
                    }
                }
            }
            if next_base.is_some() {
                let next_base_: &AlignmentStructureBase = next_base.unwrap();
                if next_base_.get_reference_chromosome_id().unwrap() == last_base.get_reference_chromosome_id().unwrap() &&
                    *next_base_.get_reference_strand().as_ref().unwrap() == *last_base.get_reference_strand().as_ref().unwrap() &&
                    next_base_.get_reference_position().unwrap().abs_diff(last_base.get_reference_position().unwrap()) <= 1 &&
                    next_base_.get_context().as_ref().unwrap().clone() == AlignmentStructureBaseContext::Exonic {
                    if last_base.get_context().as_ref().unwrap() == &AlignmentStructureBaseContext::Intronic {
                        variant_type = VariantType::IntronRetention;
                    }
                    if last_base.get_context().as_ref().unwrap() == &AlignmentStructureBaseContext::Intergenic {
                        variant_type = VariantType::UTRExtension;
                    }
                }
            }

            let graph_operation: GraphOperation = GraphOperation::new(
                first_base.get_reference_chromosome_id().unwrap(),
                first_base.get_reference_position().unwrap(),
                first_base.get_reference_strand().as_ref().unwrap().clone(),
                GraphOperationType::Include,
                last_base.get_reference_chromosome_id().unwrap(),
                last_base.get_reference_position().unwrap(),
                last_base.get_reference_strand().as_ref().unwrap().clone(),
                GraphOperationType::Include,
                "".into(),
                variant_type
            );

            variant_records.push(
                VariantRecord::new(
                    self.read_id,
                    first_base.get_read_position(),
                    last_base.get_read_position(),
                    graph_operation
                )
            );
        }

        variant_records
    }

    /// Identifies variant records in self.events based on the AlignmentStructureEventKind.
    ///
    /// This function identifies the following variant types:
    /// - Breakpoint
    /// - Deletion
    /// - Translocation
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality.
    /// * `min_base_quality`: Minimum base quality.
    ///
    /// # Returns
    /// * Vector of VariantRecord objects.
    fn identify_event_kind_variant_records(&self, min_mapping_quality: u32) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        for ((read_position_1, read_position_2), event) in self.get_events().iter() {
            let base_1: &AlignmentStructureBase = self.get_base(*read_position_1);
            let base_2: &AlignmentStructureBase = self.get_base(*read_position_2);

            // Check if the minimum mapping quality is met
            if base_1.get_mapping_quality().unwrap() < min_mapping_quality ||
                base_2.get_mapping_quality().unwrap() < min_mapping_quality {
                continue;
            }

            match event.get_kind() {
                AlignmentStructureEventKind::Breakpoint => {
                    if event.get_context().is_some() {
                        // Let the event context decide the variant type
                        continue;
                    }

                    // Fetch any insertion sequence between the breakpoints
                    let mut insertion_bases: Vec<&AlignmentStructureBase> = Vec::new();
                    let mut insertion_sequence: String = "".to_string();
                    if read_position_1.abs_diff(*read_position_2) > 1 {
                        for i in read_position_1 + 1..=read_position_2 - 1 {
                            insertion_bases.push(self.get_base(i));
                        }
                        insertion_sequence = insertion_bases
                            .iter()
                            .map(|b| b.get_nucleotide().as_str())
                            .collect::<Vec<_>>()
                            .join("");
                        let insertion_is_embedded_insertion: HashSet<bool> = insertion_bases
                            .iter()
                            .map(|b| b.is_embedded_insertion())
                            .collect();
                        assert!(
                            insertion_is_embedded_insertion.len() == 1 &&
                            insertion_is_embedded_insertion.contains(&true),
                            "All insertion nucleotides are expected to be labeled as embedded insertions."
                        );
                    }

                    let variant_type: VariantType = if base_1.get_reference_chromosome_id().unwrap() == base_2.get_reference_chromosome_id().unwrap() {
                        VariantType::Breakpoint
                    } else {
                        VariantType::Translocation
                    };

                    let graph_operation: GraphOperation = if
                        event.get_prev_read_position() == base_1.get_read_position() &&
                        event.get_next_read_position() == base_2.get_read_position() {
                        GraphOperation::new(
                            base_1.get_reference_chromosome_id().unwrap(),
                            base_1.get_reference_position().unwrap(),
                            base_1.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_prev_graph_operation_type().clone(),
                            base_2.get_reference_chromosome_id().unwrap(),
                            base_2.get_reference_position().unwrap(),
                            base_2.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_next_graph_operation_type().clone(),
                            insertion_sequence.into(),
                            variant_type
                        )
                    } else if
                        event.get_prev_read_position() == base_2.get_read_position() &&
                        event.get_next_read_position() == base_1.get_read_position() {
                        GraphOperation::new(
                            base_2.get_reference_chromosome_id().unwrap(),
                            base_2.get_reference_position().unwrap(),
                            base_2.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_prev_graph_operation_type().clone(),
                            base_1.get_reference_chromosome_id().unwrap(),
                            base_1.get_reference_position().unwrap(),
                            base_1.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_next_graph_operation_type().clone(),
                            insertion_sequence.into(),
                            variant_type
                        )
                    } else {
                        panic!("Mismatch between alignment event and bases.");
                    };

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            *read_position_1,
                            *read_position_2,
                            graph_operation
                        )
                    );
                },
                AlignmentStructureEventKind::Deletion => {
                    let graph_operation: GraphOperation = if
                        event.get_prev_read_position() == base_1.get_read_position() &&
                        event.get_next_read_position() == base_2.get_read_position() {
                        GraphOperation::new(
                            base_1.get_reference_chromosome_id().unwrap(),
                            base_1.get_reference_position().unwrap(),
                            base_1.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_prev_graph_operation_type().clone(),
                            base_2.get_reference_chromosome_id().unwrap(),
                            base_2.get_reference_position().unwrap(),
                            base_2.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_next_graph_operation_type().clone(),
                            "".into(),
                            VariantType::Deletion
                        )
                    } else if
                        event.get_prev_read_position() == base_2.get_read_position() &&
                        event.get_next_read_position() == base_1.get_read_position() {
                        GraphOperation::new(
                            base_2.get_reference_chromosome_id().unwrap(),
                            base_2.get_reference_position().unwrap(),
                            base_2.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_prev_graph_operation_type().clone(),
                            base_1.get_reference_chromosome_id().unwrap(),
                            base_1.get_reference_position().unwrap(),
                            base_1.get_reference_strand().as_ref().unwrap().clone(),
                            event.get_next_graph_operation_type().clone(),
                            "".into(),
                            VariantType::Deletion
                        )
                    } else {
                        panic!("Mismatch between alignment event and bases.");
                    };

                    variant_records.push(
                        VariantRecord::new(
                            self.read_id,
                            *read_position_1,
                            *read_position_2,
                            graph_operation
                        )
                    );
                },
                _ => {
                    // Do nothing
                }
            }
        }
        variant_records
    }

    /// Identifies variant records in self.events based on the AlignmentStructureEventContext.
    ///
    /// This function identifies the following variant types:
    /// - Circular RNA
    /// - Exon truncation
    /// - Fusion gene
    ///
    /// # Arguments
    /// * `min_mapping_quality`: Minimum mapping quality.
    /// * `min_base_quality`: Minimum base quality.
    ///
    /// # Returns
    /// * Vector of VariantRecord objects.
    fn identify_event_context_variant_records(
        &self,
        min_mapping_quality: u32,
        gene_annotator: &(impl GeneAnnotator + Sync)
    ) -> Vec<VariantRecord> {
        // Step 1. Record variant records
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        let last_base_position: u32 = self.get_bases_length() - 1;
        for ((read_position_1, read_position_2), event) in self.get_events().iter() {
            let (base_1, base_2) = if *read_position_1 == 0 && *read_position_2 == 0 {
                (self.get_base(*read_position_2), self.get_base(*read_position_2))
            } else if *read_position_1 == last_base_position && *read_position_2 == last_base_position {
                (self.get_base(*read_position_1), self.get_base(*read_position_1))
            } else {
                (self.get_base(*read_position_1), self.get_base(*read_position_2))
            };

            // Check if the minimum mapping quality is met
            if base_1.get_mapping_quality().unwrap() < min_mapping_quality ||
                base_2.get_mapping_quality().unwrap() < min_mapping_quality {
                continue;
            }

            if event.get_context().is_some() {
                match event.get_context().unwrap() {
                    AlignmentStructureEventContext::BackSplicing => {
                        // Fetch any insertion sequence between the breakpoints
                        let mut insertion_alignment_bases: Vec<&AlignmentStructureBase> = Vec::new();
                        let mut insertion_sequence: String = "".to_string();
                        if read_position_1.abs_diff(*read_position_2) > 1 {
                            for i in read_position_1 + 1..=read_position_2 - 1 {
                                insertion_alignment_bases.push(self.get_base(i));
                            }
                            insertion_sequence = insertion_alignment_bases
                                .iter()
                                .map(|b| b.get_nucleotide().as_str())
                                .collect::<Vec<_>>()
                                .join("");
                            let insertion_is_embedded_insertion: HashSet<bool> = insertion_alignment_bases
                                .iter()
                                .map(|b| b.is_embedded_insertion())
                                .collect();
                            assert!(
                                insertion_is_embedded_insertion.len() == 1 && insertion_is_embedded_insertion.contains(&true),
                                "All insertion nucleotides are expected to be labeled as embedded insertions."
                            );
                        }

                        let graph_operation: GraphOperation = if
                            event.get_prev_read_position() == base_1.get_read_position() &&
                            event.get_next_read_position() == base_2.get_read_position() {
                            GraphOperation::new(
                                base_1.get_reference_chromosome_id().unwrap(),
                                base_1.get_reference_position().unwrap(),
                                base_1.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_prev_graph_operation_type().clone(),
                                base_2.get_reference_chromosome_id().unwrap(),
                                base_2.get_reference_position().unwrap(),
                                base_2.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_next_graph_operation_type().clone(),
                                insertion_sequence.into(),
                                VariantType::CircularRNA
                            )
                        } else if
                            event.get_prev_read_position() == base_2.get_read_position() &&
                            event.get_next_read_position() == base_1.get_read_position() {
                            GraphOperation::new(
                                base_2.get_reference_chromosome_id().unwrap(),
                                base_2.get_reference_position().unwrap(),
                                base_2.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_prev_graph_operation_type().clone(),
                                base_1.get_reference_chromosome_id().unwrap(),
                                base_1.get_reference_position().unwrap(),
                                base_1.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_next_graph_operation_type().clone(),
                                insertion_sequence.into(),
                                VariantType::CircularRNA
                            )
                        } else {
                            panic!("Mismatch between alignment event and bases.");
                        };

                        variant_records.push(
                            VariantRecord::new(
                                self.read_id,
                                *read_position_1,
                                *read_position_2,
                                graph_operation
                            )
                        );
                    },
                    AlignmentStructureEventContext::FusionGene => {
                        // Fetch any insertion sequence between the breakpoints
                        let mut insertion_alignment_bases: Vec<&AlignmentStructureBase> = Vec::new();
                        let mut insertion_sequence: String = "".to_string();
                        if read_position_1.abs_diff(*read_position_2) > 1 {
                            for i in read_position_1 + 1..=read_position_2 - 1 {
                                insertion_alignment_bases.push(self.get_base(i));
                            }
                            insertion_sequence = insertion_alignment_bases
                                .iter()
                                .map(|b| b.get_nucleotide().as_str())
                                .collect::<Vec<_>>()
                                .join("");
                            let insertion_is_embedded_insertion: HashSet<bool> = insertion_alignment_bases
                                .iter()
                                .map(|b| b.is_embedded_insertion())
                                .collect();
                            assert!(
                                insertion_is_embedded_insertion.len() == 1 && insertion_is_embedded_insertion.contains(&true),
                                "All insertion nucleotides are expected to be labeled as embedded insertions."
                            );
                        }

                        let graph_operation: GraphOperation = if
                            event.get_prev_read_position() == base_1.get_read_position() &&
                            event.get_next_read_position() == base_2.get_read_position() {
                            GraphOperation::new(
                                base_1.get_reference_chromosome_id().unwrap(),
                                base_1.get_reference_position().unwrap(),
                                base_1.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_prev_graph_operation_type().clone(),
                                base_2.get_reference_chromosome_id().unwrap(),
                                base_2.get_reference_position().unwrap(),
                                base_2.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_next_graph_operation_type().clone(),
                                insertion_sequence.into(),
                                VariantType::FusionGene
                            )
                        } else if
                            event.get_prev_read_position() == base_2.get_read_position() &&
                            event.get_next_read_position() == base_1.get_read_position() {
                            GraphOperation::new(
                                base_2.get_reference_chromosome_id().unwrap(),
                                base_2.get_reference_position().unwrap(),
                                base_2.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_prev_graph_operation_type().clone(),
                                base_1.get_reference_chromosome_id().unwrap(),
                                base_1.get_reference_position().unwrap(),
                                base_1.get_reference_strand().as_ref().unwrap().clone(),
                                event.get_next_graph_operation_type().clone(),
                                insertion_sequence.into(),
                                VariantType::FusionGene
                            )
                        } else {
                            panic!("Mismatch between alignment event and bases.");
                        };

                        variant_records.push(
                            VariantRecord::new(
                                self.read_id,
                                *read_position_1,
                                *read_position_2,
                                graph_operation
                            )
                        );

                        // Record exon truncations
                        for reference_bases in event.get_skipped_reference_bases_clusters() {
                            // Identify start and end positions of exon truncation
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
                                    *read_position_1,
                                    *read_position_2,
                                    graph_operation
                                )
                            );
                        }
                    },
                    AlignmentStructureEventContext::NonCanonicalSplicing => {
                        // Record exon truncations
                        for reference_bases in event.get_skipped_reference_bases_clusters() {
                            // Identify start and end positions of exon truncation
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
                                    *read_position_1,
                                    *read_position_2,
                                    graph_operation
                                )
                            );
                        }
                    },
                    _ => {
                        // Do nothing
                    }
                }
            }
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
