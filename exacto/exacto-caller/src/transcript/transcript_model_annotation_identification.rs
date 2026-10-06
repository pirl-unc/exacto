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
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;
use crate::transcript::transcript_model_annotation::SkippedReferenceRun;


struct ReferenceTranscriptExons<'a> {
    transcript: &'a Transcript,
    chromosome_id: ReferenceChromosomeID,

    /// In transcript order: `Transcript::get_sorted_exons`.
    exons: Vec<&'a Exon>,

    /// Indices into `exons`, by first base. Empty where two exons of the transcript overlap.
    exons_by_start: Vec<usize>,

    /// First and last base of the exons: `u32::MAX` and 0 where they hold no base.
    start: ReferencePosition,
    end: ReferencePosition
}

impl<'a> ReferenceTranscriptExons<'a> {
    fn new(transcript: &'a Transcript, chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>) -> Self {
        let chromosome_id: ReferenceChromosomeID = *chromosome_names_map.get_by_left(&*transcript.chromosome).unwrap();
        let exons: Vec<&Exon> = transcript.get_sorted_exons();
        for exon in exons.iter() {
            assert_eq!(exon.strand, exons[0].strand, "Reference strand mismatch.");
        }
        let mut exons_by_start: Vec<usize> = (0..exons.len()).collect();
        exons_by_start.sort_by_key(|&index| exons[index].start);
        let is_disjoint: bool = exons_by_start
            .windows(2)
            .all(|pair| exons[pair[0]].end < exons[pair[1]].start);
        if !is_disjoint {
            exons_by_start.clear();
        }
        let start: ReferencePosition = exons.iter().filter(|exon| exon.start <= exon.end).map(|exon| exon.start).min().unwrap_or(u32::MAX);
        let end: ReferencePosition = exons.iter().filter(|exon| exon.start <= exon.end).map(|exon| exon.end).max().unwrap_or(0);

        Self {
            transcript,
            chromosome_id,
            exons,
            exons_by_start,
            start,
            end
        }
    }

    /// The index of the exon holding `position`: the last one in transcript order where two
    /// overlap.
    fn find_exon(&self, position: ReferencePosition) -> Option<usize> {
        if self.exons_by_start.is_empty() {
            return self.exons.iter().rposition(|exon| exon.start <= position && position <= exon.end);
        }
        let i: usize = self.exons_by_start.partition_point(|&index| self.exons[index].start <= position);
        if i == 0 {
            return None;
        }
        let index: usize = self.exons_by_start[i - 1];
        (position <= self.exons[index].end).then_some(index)
    }

    /// The index of the exon holding the base of the reference at `position`.
    fn find_exon_at(&self, chromosome_id: ReferenceChromosomeID, position: ReferencePosition, strand: &Strand) -> Option<usize> {
        if self.chromosome_id != chromosome_id {
            return None;
        }
        self.find_exon(position).filter(|&index| self.exons[index].strand == *strand)
    }

    /// The ID of the exon a base of the read at `position` is annotated with.
    fn get_exon_id(&self, position: ReferencePosition) -> Option<ReferenceExonID> {
        if !self.exons_by_start.is_empty() {
            return self.find_exon(position).map(|index| self.exons[index].exon_id.clone());
        }
        self.transcript
            .exons
            .values()
            .find(|exon| overlaps(position as isize, position as isize, exon.start as isize, exon.end as isize))
            .map(|exon| exon.exon_id.clone())
    }

    /// The introns between the exons, as (first base, last base).
    fn get_introns(&self) -> Vec<(ReferencePosition, ReferencePosition)> {
        let is_forward: bool = self.transcript.strand == Strand::Forward;
        let exons: Vec<&&Exon> = self.exons.iter().filter(|exon| exon.start <= exon.end).collect();
        let mut introns: Vec<(ReferencePosition, ReferencePosition)> = Vec::new();
        for pair in exons.windows(2) {
            // The last base of one exon and the first base of the next, walking the transcript.
            let (current, next): (ReferencePosition, ReferencePosition) = if is_forward {
                (pair[0].end, pair[1].start)
            } else {
                (pair[0].start, pair[1].end)
            };
            if current.abs_diff(next) > 1 {
                if pair[0].strand == Strand::Forward {
                    introns.push((current + 1, next - 1));
                } else {
                    introns.push((next + 1, current - 1));
                }
            }
        }
        introns
    }
}


/// The annotation of a base of the read, before its IDs are spelled out.
#[derive(Clone, Copy, PartialEq)]
enum BaseContext {
    /// (reference transcript, exon)
    Exonic(usize, usize),

    /// (reference transcript)
    Intronic(usize),

    Intergenic
}


pub fn identify_transcript_model_annotation(
    alignment_model: &AlignmentModel,
    reference_transcript_matches: &Vec<ReferenceTranscriptMatch>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap
) -> TranscriptModelAnnotation {
    // Step 1. Get the exons of the reference transcripts.
    let references: Vec<ReferenceTranscriptExons> = reference_transcript_matches
        .iter()
        .map(|reference_transcript_match| {
            let reference_transcript: &Transcript = gene_annotator
                .get_transcript(&*reference_transcript_match.get_reference_transcript_id())
                .unwrap();
            ReferenceTranscriptExons::new(reference_transcript, chromosome_names_map)
        })
        .collect();

    // Step 2. Each reference transcript should have at least 1 base that overlaps a model base.
    for reference in references.iter() {
        let has_overlap: bool = alignment_model.get_bases().iter()
            .filter(|base| matches!(
                base.get_kind(),
                AlignmentModelBaseKind::Match |
                AlignmentModelBaseKind::Mismatch |
                AlignmentModelBaseKind::Insertion
            ))
            .filter_map(|base| base.get_placement().get_coordinate())
            .any(|(reference_chromosome_id, reference_position, reference_strand)| {
                reference.find_exon_at(reference_chromosome_id, reference_position, reference_strand).is_some()
            });
        if !has_overlap {
            panic!(
                "None of the ReferenceTranscriptSequence bases for {} overlaps with any of the AlignmentModel bases.",
                reference.transcript.transcript_id
            );
        }
    }

    // Step 3. Identify the context of each AlignmentModelBase. A base of the reference that
    // two reference transcripts hold belongs to the earlier one: the matches come best
    // first. The aligned bases come in read order, so the annotation is appended run by run.
    let mut annotation: TranscriptModelAnnotation = TranscriptModelAnnotation::new();
    // Vec<(read position, reference position)> by reference transcript: 
    // the first and the last base of the read annotated with it, in reference order.
    let mut first_bases: Vec<Option<(ReadPosition, ReferencePosition)>> = vec![None; references.len()];
    let mut last_bases: Vec<Option<(ReadPosition, ReferencePosition)>> = vec![None; references.len()];
    // (first read position, last read position, context)
    let mut run: Option<(ReadPosition, ReadPosition, BaseContext)> = None;
    let push_run = |annotation: &mut TranscriptModelAnnotation, (read_start, read_end, context): (ReadPosition, ReadPosition, BaseContext)| {
        let base_annotation: TranscriptModelBaseAnnotation = match context {
            BaseContext::Exonic(reference, _) => {
                let (_, reference_position, _) = alignment_model.get_base(read_start).get_placement().get_coordinate().unwrap();
                let reference_exon_id: Option<ReferenceExonID> = references[reference].get_exon_id(reference_position);
                assert_eq!(reference_exon_id.is_some(), true);
                TranscriptModelBaseAnnotation {
                    context: AlignmentModelBaseContext::Exonic,
                    reference_gene_id: Some(references[reference].transcript.gene_id.clone()),
                    reference_transcript_id: Some(references[reference].transcript.transcript_id.clone()),
                    reference_exon_id: reference_exon_id
                }
            },
            BaseContext::Intronic(reference) => TranscriptModelBaseAnnotation {
                context: AlignmentModelBaseContext::Intronic,
                reference_gene_id: Some(references[reference].transcript.gene_id.clone()),
                reference_transcript_id: Some(references[reference].transcript.transcript_id.clone()),
                reference_exon_id: None
            },
            BaseContext::Intergenic => TranscriptModelBaseAnnotation {
                context: AlignmentModelBaseContext::Intergenic,
                reference_gene_id: None,
                reference_transcript_id: None,
                reference_exon_id: None
            }
        };
        annotation.push_bases(read_start, read_end, base_annotation);
    };
    for i in 0..alignment_model.num_bases() {
        let base: &AlignmentModelBase = alignment_model.get_base(i);
        if !matches!(
            base.get_kind(),
            | AlignmentModelBaseKind::Match
            | AlignmentModelBaseKind::Mismatch
            | AlignmentModelBaseKind::Insertion
        ) {
            continue;
        }
        assert!(base.get_placement().is_placed() == true, "Aligned base at read position {} is unplaced.", i);
        let (chromosome_id, position, strand) = base.get_placement().get_coordinate().unwrap();
        let exonic: Option<(usize, usize)> = references
            .iter()
            .enumerate()
            .find_map(|(index, reference)| reference.find_exon_at(chromosome_id, position, strand).map(|exon| (index, exon)));
        let context: BaseContext = match exonic {
            Some((reference, exon)) if !references[reference].exons_by_start.is_empty() => BaseContext::Exonic(reference, exon),
            // Where exons of the transcript overlap, every base is its own run.
            Some((reference, _)) => BaseContext::Exonic(reference, usize::MAX),
            // Intronic: inside the span of a transcript on the same chromosome and strand, as
            // the exon test above asks.
            None => references
                .iter()
                .position(|reference| {
                    reference.chromosome_id == chromosome_id
                        && reference.transcript.strand == *strand
                        && position >= reference.start
                        && position <= reference.end
                })
                .map_or(BaseContext::Intergenic, BaseContext::Intronic)
        };
        if let BaseContext::Exonic(reference, _) | BaseContext::Intronic(reference) = context {
            if first_bases[reference].map_or(true, |(_, first)| position < first) {
                first_bases[reference] = Some((base.get_read_position(), position));
            }
            if last_bases[reference].map_or(true, |(_, last)| position >= last) {
                last_bases[reference] = Some((base.get_read_position(), position));
            }
        }

        let read_position: ReadPosition = base.get_read_position();
        run = match run {
            Some((read_start, read_end, run_context))
                if read_end + 1 == read_position
                    && run_context == context
                    && !matches!(context, BaseContext::Exonic(_, usize::MAX)) =>
            {
                Some((read_start, read_position, run_context))
            },
            Some(finished) => {
                push_run(&mut annotation, finished);
                Some((read_position, read_position, context))
            },
            None => Some((read_position, read_position, context))
        };
    }
    if let Some(finished) = run {
        push_run(&mut annotation, finished);
    }

    // Step 4. Identify reference transcript introns.
    // HashMap<chromosome ID, HashSet<(start, end)>>
    let mut reference_transcripts_introns_map: HashMap<ReferenceChromosomeID, HashSet<(ReferencePosition, ReferencePosition)>> = HashMap::new();
    for reference in references.iter() {
        for (start, end) in reference.get_introns() {
            reference_transcripts_introns_map
                .entry(reference.chromosome_id)
                .or_insert_with(HashSet::new)
                .insert((start, end));
        }
    }

    // Step 5. Identify the context of each AlignmentModelEvent.
    let event_keys: Vec<(ReadPosition, ReadPosition)> = alignment_model.get_events().keys().cloned().collect();
    for (read_position_1, read_position_2) in event_keys {
        let base_1: &AlignmentModelBase = alignment_model.get_base(read_position_1);
        let base_2: &AlignmentModelBase = alignment_model.get_base(read_position_2);
        let alignment_event: &AlignmentModelEvent = alignment_model.get_event(read_position_1, read_position_2).unwrap();
        let (reference_chromosome_id_1, reference_position_1, reference_strand_1) = base_1.get_placement().get_coordinate().unwrap();
        let (reference_chromosome_id_2, reference_position_2, _) = base_2.get_placement().get_coordinate().unwrap();
        let chromosome_name_1: &ReferenceChromosomeName = chromosome_names_map.get_by_right(&reference_chromosome_id_1).unwrap();
        let chromosome_name_2: &ReferenceChromosomeName = chromosome_names_map.get_by_right(&reference_chromosome_id_2).unwrap();
        // Each side of the event lies in a gene and the two sides share none. A junction
        // inside one gene is not a fusion because another gene overlaps a flank.
        let joins_genes = || -> bool {
            let gene_ids_1: HashSet<ReferenceGeneID> = gene_annotator.get_gene_ids_at_locus(&**chromosome_name_1, reference_position_1).into_iter().collect();
            let gene_ids_2: HashSet<ReferenceGeneID> = gene_annotator.get_gene_ids_at_locus(&**chromosome_name_2, reference_position_2).into_iter().collect();
            !gene_ids_1.is_empty() && !gene_ids_2.is_empty() && gene_ids_1.is_disjoint(&gene_ids_2)
        };
        match *alignment_event.get_kind() {
            AlignmentModelEventKind::Splicing => {
                let reference_chromosome_id: ReferenceChromosomeID = reference_chromosome_id_1;
                let reference_strand: Strand = reference_strand_1.clone();
                let mut reference_start: ReferencePosition = reference_position_1;
                let mut reference_end: ReferencePosition = reference_position_2;
                if reference_strand == Strand::Forward {
                    reference_start = reference_start + 1;
                    reference_end = reference_end - 1;
                } else {
                    reference_start = reference_start - 1;
                    reference_end = reference_end + 1;
                }
                if reference_transcripts_introns_map.contains_key(&reference_chromosome_id) {
                    let reference_transcripts_introns: &HashSet<(ReferencePosition, ReferencePosition)> = reference_transcripts_introns_map.get(&reference_chromosome_id).unwrap();
                    if reference_transcripts_introns.contains(&(reference_start, reference_end))
                        || reference_transcripts_introns.contains(&(reference_end, reference_start)) {
                        annotation.set_event_context(
                            read_position_1,
                            read_position_2,
                            AlignmentModelEventContext::CanonicalSplicing
                        );
                    } else {
                        if joins_genes() {
                            annotation.set_event_context(
                                read_position_1,
                                read_position_2,
                                AlignmentModelEventContext::FusionGene
                            );
                        } else {
                            annotation.set_event_context(
                                read_position_1,
                                read_position_2,
                                AlignmentModelEventContext::NonCanonicalSplicing
                            );
                        }
                    }
                } else {
                    annotation.set_event_context(
                        read_position_1,
                        read_position_2,
                        AlignmentModelEventContext::NonCanonicalSplicing
                    );
                }
            },
            AlignmentModelEventKind::Breakpoint => {
                if reference_chromosome_id_1 == reference_chromosome_id_2 {
                    // Build sets of (chromosome_id, reference_position, strand) tuples for
                    // bases flanking the breakpoint. Bases without a reference assignment
                    // (soft-clipped, inserted, etc.) are `BasePlacement::Unplaced`
                    // and are excluded. A read that folds back onto the other strand of the
                    // same place shares positions but not the strand, and is no back-splice.
                    let left_bases_set: HashSet<(ReferenceChromosomeID, ReferencePosition, Strand)> = (0..=read_position_1)
                        .filter_map(|i| {
                            alignment_model.get_base(i)
                                .get_placement()
                                .get_coordinate()
                                .map(|(chromosome_id, position, strand)| (chromosome_id, position, strand.clone()))
                        })
                        .collect();
                    let right_bases_set: HashSet<(ReferenceChromosomeID, ReferencePosition, Strand)> = (read_position_2..alignment_model.num_bases())
                        .filter_map(|i| {
                            alignment_model.get_base(i)
                                .get_placement()
                                .get_coordinate()
                                .map(|(chromosome_id, position, strand)| (chromosome_id, position, strand.clone()))
                        })
                        .collect();
                    if left_bases_set.is_disjoint(&right_bases_set) == false {
                        annotation.set_event_context(
                            read_position_1,
                            read_position_2,
                            AlignmentModelEventContext::BackSplicing
                        );
                    } else {
                        if joins_genes() {
                            annotation.set_event_context(
                                read_position_1,
                                read_position_2,
                                AlignmentModelEventContext::FusionGene
                            );
                        } else {
                            annotation.set_event_context(
                                read_position_1,
                                read_position_2,
                                AlignmentModelEventContext::NonCanonicalSplicing
                            );
                        }
                    }
                } else {
                    if joins_genes() {
                        annotation.set_event_context(
                            read_position_1,
                            read_position_2,
                            AlignmentModelEventContext::FusionGene
                        );
                    }
                }
            },
            _ => {
                // Do nothing.
            }
        }
    }

    // Step 6. Get the bases of the read that anchor each reference transcript: the ends of
    // the events annotated with it, and the first and the last base annotated with it.
    // Sorted by reference position.
    // Vec<Vec<(read position, reference position)>>, by reference transcript
    let mut anchors: Vec<Vec<(ReadPosition, ReferencePosition)>> = vec![Vec::new(); references.len()];
    let reference_indices: HashMap<&str, usize> = references
        .iter()
        .enumerate()
        .map(|(index, reference)| (&*reference.transcript.transcript_id, index))
        .collect();
    for (read_position_1, read_position_2) in alignment_model.get_events().keys() {
        for read_position in [*read_position_1, *read_position_2] {
            // Unannotated bases and intergenic bases (annotated, but with no
            // reference transcript ID) are not exonic boundaries.
            let Some(base_annotation) = annotation.get_base(read_position) else {
                continue;
            };
            let Some(reference_transcript_id) = base_annotation.reference_transcript_id.as_ref() else {
                continue;
            };
            // Only placed bases are annotated.
            let (_, reference_position, _) = alignment_model
                .get_base(read_position)
                .get_placement()
                .get_coordinate()
                .unwrap();
            anchors[reference_indices[&**reference_transcript_id]].push((read_position, reference_position));
        }
    }
    for (index, reference) in references.iter().enumerate() {
        match (first_bases[index], last_bases[index]) {
            (Some(first), Some(last)) => {
                anchors[index].push(first);
                anchors[index].push(last);
            },
            _ => {
                // No base of the read is annotated with the reference transcript: its
                // anchors are the first and the last base of the read on its exons. Step 2
                // guarantees at least one.
                assert!(anchors[index].is_empty());
                let mut first: Option<(ReadPosition, ReferencePosition)> = None;
                let mut last: Option<(ReadPosition, ReferencePosition)> = None;
                for base in alignment_model.get_bases().iter().filter(|base| matches!(
                    base.get_kind(),
                    AlignmentModelBaseKind::Match |
                    AlignmentModelBaseKind::Mismatch |
                    AlignmentModelBaseKind::Insertion
                )) {
                    let Some((chromosome_id, position, strand)) = base.get_placement().get_coordinate() else {
                        continue;
                    };
                    if reference.find_exon_at(chromosome_id, position, strand).is_none() {
                        continue;
                    }
                    if first.map_or(true, |(_, reference_position)| position < reference_position) {
                        first = Some((base.get_read_position(), position));
                    }
                    if last.map_or(true, |(_, reference_position)| position >= reference_position) {
                        last = Some((base.get_read_position(), position));
                    }
                }
                anchors[index].push(first.unwrap());
                anchors[index].push(last.unwrap());
            }
        }
    }
    for reference_anchors in anchors.iter_mut() {
        reference_anchors.sort_unstable_by_key(|&(read_position, reference_position)| (reference_position, read_position));
        reference_anchors.dedup();
    }

    // Step 7. Identify the bases of the reference transcripts the read skipped. They are the
    // exons minus the placed bases of the read, run by run, and each run goes to the event
    // of the anchor closest to it.
    // BTreeMap<(chromosome ID, strand), Vec<(start, end, reference transcript, exon)>>
    let mut exons: BTreeMap<(ReferenceChromosomeID, &Strand), Vec<(ReferencePosition, ReferencePosition, usize, usize)>> = BTreeMap::new();
    for (index, reference) in references.iter().enumerate() {
        for (exon_index, exon) in reference.exons.iter().enumerate() {
            if exon.start <= exon.end {
                exons
                    .entry((reference.chromosome_id, &exon.strand))
                    .or_default()
                    .push((exon.start, exon.end, index, exon_index));
            }
        }
    }
    // HashMap<(chromosome ID, strand), Vec<reference position>>
    let mut placed_positions: HashMap<(ReferenceChromosomeID, &Strand), Vec<ReferencePosition>> = HashMap::new();
    for base in alignment_model.get_bases() {
        if matches!(
            base.get_kind(),
            AlignmentModelBaseKind::Match
            | AlignmentModelBaseKind::Mismatch
            | AlignmentModelBaseKind::Insertion
        ) {
            let (chromosome_id, position, strand) = base.get_placement().get_coordinate().unwrap();
            placed_positions.entry((chromosome_id, strand)).or_default().push(position);
        }
    }
    for positions in placed_positions.values_mut() {
        positions.sort_unstable();
        positions.dedup();
    }

    // The other end of the anchor's event that a skipped base goes to. A base can end two events,
    // a deletion on one side of it and an intron on the other, and `get_event_at` holds only one:
    // the skipped base goes to the event whose other end lies on its side of the anchor, and to
    // the one `get_event_at` holds when none does.
    let event_beside = |
        anchor_read_position: ReadPosition,
        anchor_reference_position: ReferencePosition,
        chromosome_id: ReferenceChromosomeID,
        reference_position: ReferencePosition
    | -> Option<ReadPosition> {
        let indexed: Option<ReadPosition> = alignment_model.get_event_at(anchor_read_position).map(|event| {
            if event.get_prev_read_position() == anchor_read_position {
                event.get_next_read_position()
            } else {
                event.get_prev_read_position()
            }
        });
        let mut beside: Vec<ReadPosition> = alignment_model
            .get_events()
            .keys()
            .filter_map(|&(read_position_1, read_position_2)| {
                if read_position_1 == anchor_read_position {
                    Some(read_position_2)
                } else if read_position_2 == anchor_read_position {
                    Some(read_position_1)
                } else {
                    None
                }
            })
            .filter(|&other_read_position| match alignment_model.get_base(other_read_position).get_placement().get_coordinate() {
                Some((other_chromosome_id, other_reference_position, _)) => {
                    other_chromosome_id == chromosome_id
                        && other_reference_position != anchor_reference_position
                        && (other_reference_position < anchor_reference_position) == (reference_position < anchor_reference_position)
                },
                None => false
            })
            .collect();
        beside.sort_unstable();
        match indexed {
            Some(other_read_position) if beside.contains(&other_read_position) => Some(other_read_position),
            _ => beside.first().copied().or(indexed)
        }
    };

    // HashMap<(read position of the anchor, whether the skipped base lies before it), (read position 1, read position 2) of its event>
    let mut anchor_events: HashMap<(ReadPosition, bool), (ReadPosition, ReadPosition)> = HashMap::new();
    for (&(chromosome_id, strand), exons) in exons.iter() {
        let no_positions: Vec<ReferencePosition> = Vec::new();
        let positions: &Vec<ReferencePosition> = placed_positions.get(&(chromosome_id, strand)).unwrap_or(&no_positions);

        // Cut the exons where one starts or stops, so every piece has one owner: the exon
        // of the earliest reference transcript holding it.
        let mut cuts: Vec<u64> = exons
            .iter()
            .flat_map(|&(start, end, _, _)| [start as u64, end as u64 + 1])
            .collect();
        cuts.sort_unstable();
        cuts.dedup();
        for cut in cuts.windows(2) {
            let (piece_start, piece_end): (ReferencePosition, ReferencePosition) = (cut[0] as ReferencePosition, (cut[1] - 1) as ReferencePosition);
            let Some(&(_, _, index, exon_index)) = exons
                .iter()
                .filter(|&&(start, end, _, _)| start <= piece_start && piece_end <= end)
                .min_by_key(|&&(_, _, index, exon_index)| (index, Reverse(exon_index))) else {
                continue;
            };
            let reference: &ReferenceTranscriptExons = &references[index];
            let exon: &Exon = reference.exons[exon_index];
            let is_forward: bool = reference.transcript.strand == Strand::Forward;

            // The stretches of the piece between the placed bases of the read.
            let first: usize = positions.partition_point(|&position| position < piece_start);
            let last: usize = positions.partition_point(|&position| position <= piece_end);
            let mut skipped: Vec<(ReferencePosition, ReferencePosition)> = Vec::new();
            let mut next_start: u64 = piece_start as u64;
            for &position in positions[first..last].iter() {
                if next_start < position as u64 {
                    skipped.push((next_start as ReferencePosition, position - 1));
                }
                next_start = position as u64 + 1;
            }
            if next_start <= piece_end as u64 {
                skipped.push((next_start as ReferencePosition, piece_end));
            }

            for (skipped_start, skipped_end) in skipped {
                let sequence: &[u8] = fasta_map
                    .get_sequence(&*reference.transcript.chromosome, skipped_start as usize, skipped_end as usize)
                    .as_bytes();
                // (event, first base, nucleotides)
                let mut run: Option<((ReadPosition, ReadPosition), ReferencePosition, String)> = None;
                for reference_position in skipped_start..=skipped_end {
                    // Identify the closest base.
                    let vec: &Vec<(ReadPosition, ReferencePosition)> = &anchors[index];
                    let (closest_read_position, closest_reference_position): (ReadPosition, ReferencePosition) = match vec.binary_search_by_key(&reference_position, |&(_, base_reference_position)| base_reference_position) {
                        Ok(idx) => {
                            vec[idx]
                        },
                        Err(idx) => {
                            if idx == 0 {
                                vec[0]
                            } else if idx == vec.len() {
                                vec[vec.len() - 1]
                            } else {
                                let reference_position_1: ReferencePosition = vec[idx-1].1;
                                let reference_position_2: ReferencePosition = vec[idx].1;
                                if reference_position.abs_diff(reference_position_1) <= reference_position.abs_diff(reference_position_2) {
                                    vec[idx-1]
                                } else {
                                    vec[idx]
                                }
                            }
                        }
                    };

                    // Identify the closest event.
                    let is_before: bool = reference_position < closest_reference_position;
                    let event: (ReadPosition, ReadPosition) = *anchor_events.entry((closest_read_position, is_before)).or_insert_with(|| {
                        if closest_read_position == 0 || closest_read_position == alignment_model.num_bases() - 1 {
                            // The skipped reference bases lie beyond the boundaries of the alignment structure
                            // for the given reference transcript.
                            annotation.set_boundary_event(closest_read_position);
                            (closest_read_position, closest_read_position)
                        } else if let Some(closest_read_position_2) = event_beside(
                            closest_read_position,
                            closest_reference_position,
                            chromosome_id,
                            reference_position
                        ) {
                            assert!(
                                alignment_model.get_event(closest_read_position, closest_read_position_2).is_some(),
                                "Event does not exist between bases {} and {}", closest_read_position, closest_read_position_2
                            );
                            (closest_read_position, closest_read_position_2)
                        } else {
                            // No event exists at this position (e.g. no splicing in this read).
                            // Create a boundary event to anchor the skipped reference base.
                            annotation.set_boundary_event(closest_read_position);
                            (closest_read_position, closest_read_position)
                        }
                    });

                    // The nucleotide on the strand of the transcript.
                    let offset: usize = (reference_position - skipped_start) as usize;
                    let nucleotide_sequence: &str = std::str::from_utf8(&sequence[offset..offset + 1])
                        .expect("Failed to convert sequence to UTF-8");
                    let mut nucleotide: Nucleotide = Nucleotide::from_str(nucleotide_sequence)
                        .unwrap_or_else(|_| if nucleotide_sequence.chars().next().is_some_and(|c| c.is_ascii_lowercase()) {
                            Nucleotide::n
                        } else {
                            Nucleotide::N
                        });
                    if !is_forward {
                        nucleotide = nucleotide.complement();
                    }

                    run = match run {
                        Some((run_event, run_start, mut nucleotides)) if run_event == event => {
                            nucleotides.push_str(nucleotide.as_str());
                            Some((run_event, run_start, nucleotides))
                        },
                        finished => {
                            if let Some(((read_position_1, read_position_2), run_start, nucleotides)) = finished {
                                annotation.add_skipped_reference_run(read_position_1, read_position_2, SkippedReferenceRun {
                                    reference_chromosome_id: chromosome_id,
                                    reference_start: run_start,
                                    reference_end: reference_position - 1,
                                    reference_strand: exon.strand.clone(),
                                    reference_gene_id: Some(exon.gene_id.clone()),
                                    reference_transcript_id: Some(exon.transcript_id.clone()),
                                    reference_exon_id: Some(exon.exon_id.clone()),
                                    sequence: nucleotides.into()
                                });
                            }
                            Some((event, reference_position, nucleotide.as_str().to_string()))
                        }
                    };
                }
                if let Some(((read_position_1, read_position_2), run_start, nucleotides)) = run {
                    annotation.add_skipped_reference_run(read_position_1, read_position_2, SkippedReferenceRun {
                        reference_chromosome_id: chromosome_id,
                        reference_start: run_start,
                        reference_end: skipped_end,
                        reference_strand: exon.strand.clone(),
                        reference_gene_id: Some(exon.gene_id.clone()),
                        reference_transcript_id: Some(exon.transcript_id.clone()),
                        reference_exon_id: Some(exon.exon_id.clone()),
                        sequence: nucleotides.into()
                    });
                }
            }
        }
    }

    annotation
}
