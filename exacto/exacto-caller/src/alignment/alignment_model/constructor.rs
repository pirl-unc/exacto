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


use super::AlignmentModel;

use bstr::ByteSlice;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use regex::Regex;
use std::collections::HashMap;
use std::str::FromStr;
use std::sync::Arc;

use crate::prelude::*;


impl AlignmentModel {
    pub fn from_read(read_id: ReadID, sequence: &str, qualities: &[BaseQuality]) -> Self {
        assert_eq!(sequence.len(), qualities.len());
        let bases: Vec<AlignmentModelBase> = sequence.chars().zip(qualities).enumerate().map(|(i, (base, &quality))| {
            AlignmentModelBase::new(
                i as ReadPosition,
                Nucleotide::from_str(&base.to_string()).unwrap(),
                quality
            )
        }).collect();
        Self {
            read_id,
            bases,
            records: Vec::new(),
            events: HashMap::new(),
            events_index: HashMap::new()
        }
    }

    pub fn new(
        read_id: ReadID,
        read_sequence: &str,
        base_quality_scores: &[BaseQuality],
        bam_records: &[Arc<bam::Record>]
    ) -> Self {
        if bam_records.is_empty() {
            panic!("bam_records is empty");
        }

        let first: &str = std::str::from_utf8(bam_records[0].name()
            .unwrap()
            .as_bytes())
            .unwrap();
        for record in bam_records.iter().skip(1) {
            let name: &str = std::str::from_utf8(record.name().unwrap().as_bytes()).unwrap();
            if name != first {
                panic!("bam_records do not all belong to one read: {first} vs {name}.");
            }
        }

        // A read stored without base qualities (SAM `*`) is given Q60 on every base.
        let default_base_quality_scores: Vec<BaseQuality>;
        let base_quality_scores: &[BaseQuality] = if base_quality_scores.is_empty() {
            default_base_quality_scores = vec![60u8; read_sequence.chars().count()];
            &default_base_quality_scores
        } else {
            base_quality_scores
        };
        assert_eq!(
            read_sequence.chars().count(),
            base_quality_scores.len(),
            "read sequence and base quality scores must be the same length."
        );

        // Step 1. Locate each record's span within the read.
        let records: Vec<AlignmentRecord> = Self::build_alignment_records(
            read_id,
            &*read_sequence,
            bam_records
        );

        // Step 2. One base per read base, all unplaced.
        let mut model: AlignmentModel = Self::from_read(read_id, &*read_sequence, base_quality_scores);
        let bases: &mut Vec<AlignmentModelBase> = &mut model.bases;
        let mut events: HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent> = HashMap::new();
        let mut events_index: HashMap<ReadPosition, ReadPosition> = HashMap::new();

        // Step 3. Place bases and add local events from the CS tag.
        Self::apply_local_variants(
            read_id,
            bases,
            &mut events,
            &mut events_index,
            &records
        );

        // Step 4. Add breakpoint events from split alignments and terminal clips.
        Self::apply_breakpoints(
            read_id,
            bases,
            &mut events,
            &mut events_index,
            &records,
            read_sequence
        );

        // Step 5. Assemble and check invariants.
        model.records = records;
        model.events = events;
        model.events_index = events_index;

        assert_eq!(
            model.get_read_sequence().to_uppercase(),
            read_sequence.to_uppercase(),
            "model bases do not round-trip to the original read sequence (read ID: {}).",
            model.read_id
        );

        for ((p1, p2), _) in model.events.iter() {
            assert!(p1 <= p2, "event key must be ordered: ({p1}, {p2})");
            assert!(
                (*p2 as usize) < model.bases.len(),
                "event references read position {p2} beyond the read"
            );
        }

        model
    }

    /// Locate each BAM record's span within the original read sequence.
    fn build_alignment_records(
        read_id: ReadID,
        read_sequence: &str,
        records: &[Arc<bam::Record>]
    ) -> Vec<AlignmentRecord> {
        let mut alignment_records: Vec<AlignmentRecord> = Vec::new();

        // Step 1. Identify alignment records.
        for record in records {
            // Get the aligned sequence
            let aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(&record).into();

            // Identify all start positions between aligned sequence and original read sequence.
            let start_positions: Vec<ReadPosition> = find_substring_positions(
                &*read_sequence.to_uppercase(),
                &*aligned_sequence.to_uppercase()
            );
            if start_positions.is_empty() {
                panic!("Could not find aligned sequence for read ID: {read_id}");
            }

            // Get left and right soft-clipping information of the current record.
            let left_softclipping: (bool, u32) = get_clipping(&record, true);
            let right_softclipping: (bool, u32) = get_clipping(&record, false);

            // If there are multiple start positions, find where the aligned sequence
            // starts on the original read sequence.
            let mut found_alignment: bool = false;
            let mut read_start: ReadPosition = 0;
            let mut read_end: ReadPosition = 0;
            let reference_strand: Strand = get_alignment_strand(&record);
            for start_position in start_positions.iter() {
                // Check if the current start position aligns with the current alignment record.
                let end_position: ReadPosition = *start_position + aligned_sequence.len() as u32 - 1;
                let num_left_bases: u32 = *start_position;
                let num_right_bases: u32 = read_sequence.len() as u32 - end_position - 1;

                let mut aligned: bool = true;
                if reference_strand == Strand::Reverse {
                    if (left_softclipping.0 && left_softclipping.1 != num_right_bases)
                        || (right_softclipping.0 && right_softclipping.1 != num_left_bases) {
                        aligned = false;
                    }
                } else {
                    if (left_softclipping.0 && left_softclipping.1 != num_left_bases)
                        || (right_softclipping.0 && right_softclipping.1 != num_right_bases) {
                        aligned = false;
                    }
                }

                if aligned {
                    read_start = *start_position;
                    read_end = read_start + (aligned_sequence.len() as u32) - 1;
                    found_alignment = true;
                    break;
                }
            }

            if !found_alignment {
                panic!("Could not find aligned sequence for read ID: {read_id}.");
            }

            assert_eq!(
                &*aligned_sequence,
                read_sequence[(read_start as usize)..(read_end as usize + 1)].to_string(),
                "Aligned sequence does not match the identified part of the original read sequence."
            );

            let alignment_record: AlignmentRecord = AlignmentRecord::new(
                read_start,
                read_end,
                reference_strand,
                record.clone()
            );

            alignment_records.push(alignment_record);
        }

        // Step 2. Sort alignment records by read start position.
        alignment_records.sort_by_key(|alignment| alignment.read_start);

        alignment_records
    }

    /// Place bases and add deletion/splicing events by walking each record's CS tag.
    fn apply_local_variants(
        read_id: ReadID,
        bases: &mut Vec<AlignmentModelBase>,
        events: &mut HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,
        events_index: &mut HashMap<ReadPosition, ReadPosition>,
        alignment_records: &[AlignmentRecord]
    ) {
        // Step 1. Update alignment sequence.
        for alignment_record in alignment_records.iter() {
            Self::purge_events(
                events,
                events_index,
                alignment_record.read_start,
                alignment_record.read_end
            );

            let mut reference_position: isize = get_alignment_start_position(&alignment_record.record) as isize - 1;
            let reference_chromosome_id: ReferenceChromosomeID = alignment_record.record.reference_sequence_id().unwrap().unwrap() as ReferenceChromosomeID;
            let reference_strand: Strand = get_alignment_strand(&alignment_record.record);
            let mapping_quality: MappingQuality = get_alignment_mapping_quality(&alignment_record.record);
            let cs_tag: String = get_tag_value(&alignment_record.record, "cs")
                .unwrap_or_else(|| {
                    panic!(
                        "Read ID {} alignment record has no cs tag. Exacto needs an aligner that emits one (e.g. minimap --cs).",
                        read_id
                    )
                })
                .to_string();
            // The long form of the tag (minimap2 --cs=long) spells a match out: `=ACGT` is `:4`.
            let cs_tag: String = if cs_tag.contains('=') {
                Regex::new(r"=([A-Za-z]+)")
                    .unwrap()
                    .replace_all(&cs_tag, |captures: &regex::Captures| format!(":{}", captures[1].len()))
                    .into_owned()
            } else {
                cs_tag
            };
            let mut read_position: isize = if reference_strand == Strand::Forward {
                alignment_record.read_start as isize - 1
            } else {
                alignment_record.read_end as isize + 1
            };

            // Identify SNVs, insertions, deletions, and splicing in the CS tag.
            let mut last_aligned_read_position: Option<ReadPosition> = None;
            let mut pending_events: Vec<(AlignmentModelEventKind, ReadPosition)> = Vec::new();
            let re: Regex = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap();
            for cap in re.captures_iter(&cs_tag) {
                let token: &str = &cap[0];
                let mut chars = token.chars();
                let cs_tag_kind: CSTagKind = CSTagKind::from_str(chars.next().unwrap().to_string().as_str()).unwrap();
                let payload: &str = chars.as_str();

                match cs_tag_kind {
                    CSTagKind::Match => {
                        let length: isize = payload.parse().unwrap();
                        for _ in 0..length {
                            read_position += if reference_strand == Strand::Forward {
                                1
                            } else {
                                -1
                            };
                            reference_position += 1;
                            Self::place_base(
                                bases,
                                read_position as ReadPosition,
                                AlignmentModelBaseKind::Match,
                                reference_chromosome_id,
                                reference_position as ReferencePosition,
                                reference_strand.clone(),
                                mapping_quality,
                                None
                            );
                            Self::flush_pending_events(
                                events,
                                events_index,
                                &mut pending_events,
                                &reference_strand,
                                read_position as ReadPosition
                            );
                            last_aligned_read_position = Some(read_position as ReadPosition);
                        }
                    },
                    CSTagKind::Mismatch => {
                        let alleles: Vec<char> = payload.chars().collect();
                        assert_eq!(alleles.len(), 2, "1 reference allele and 1 alternate allele expected.");
                        let mut reference_nucleotide: Nucleotide = Nucleotide::from_str(
                            alleles[0].to_ascii_uppercase().to_string().as_str()
                        ).unwrap_or_else(|_| {
                            panic!(
                                "Read ID {} contains unrecognized reference allele in the cs tag: {}.",
                                read_id, cs_tag
                            )
                        });
                        if reference_strand == Strand::Reverse {
                            reference_nucleotide = reference_nucleotide.complement();
                        }
                        read_position += if reference_strand == Strand::Forward { 1 } else { -1 };
                        reference_position += 1;
                        Self::place_base(
                            bases,
                            read_position as ReadPosition,
                            AlignmentModelBaseKind::Mismatch,
                            reference_chromosome_id,
                            reference_position as ReferencePosition,
                            reference_strand.clone(),
                            mapping_quality,
                            Some(reference_nucleotide)
                        );
                        Self::flush_pending_events(
                            events,
                            events_index,
                            &mut pending_events,
                            &reference_strand,
                            read_position as ReadPosition
                        );
                        last_aligned_read_position = Some(read_position as ReadPosition);
                    },
                    CSTagKind::Insertion => {
                        let insertion: String = payload.to_string();
                        let length: usize = insertion.chars().count();
                        for _ in 0..length {
                            read_position += if reference_strand == Strand::Forward { 1 } else { -1 };
                            Self::place_base(
                                bases,
                                read_position as ReadPosition,
                                AlignmentModelBaseKind::Insertion,
                                reference_chromosome_id,
                                reference_position as ReferencePosition,
                                reference_strand.clone(),
                                mapping_quality,
                                None
                            );
                        }
                    }
                    CSTagKind::Deletion => {
                        let length: usize = payload.chars().count();
                        if let Some(left_flank) = last_aligned_read_position {
                            pending_events.push((AlignmentModelEventKind::Deletion, left_flank));
                            let mut deleted: Vec<Nucleotide> = payload
                                .chars()
                                .map(|c| {
                                    Nucleotide::from_str(c.to_ascii_uppercase().to_string().as_str())
                                        .unwrap_or_else(|_| {
                                            panic!(
                                                "Read ID {} contains unrecognized deleted base in the cs tag: {}.",
                                                read_id, cs_tag
                                            )
                                        })
                                })
                                .collect();
                            if reference_strand == Strand::Reverse {
                                deleted = deleted
                                    .into_iter()
                                    .rev()
                                    .map(|base| base.complement())
                                    .collect();
                            }
                            bases[left_flank as usize].set_deleted_reference_bases(deleted);
                        }
                        reference_position += length as isize;
                    }
                    CSTagKind::Splicing => {
                        let re_splicing: Regex = Regex::new(r"\d+").unwrap();
                        let caps: regex::Match<'_> = re_splicing.find(&payload).expect("No numerical value found");
                        let num_start: usize = caps.start();
                        let num_end: usize = caps.end();

                        // Extract donor splice site signal (2 letters before the number).
                        let mut donor_splice_site_signal: Box<str> = payload[num_start - 2..num_start].into();

                        // Extract acceptor splice site signal (2 letters after the number).
                        let mut acceptor_splice_site_signal: Box<str> = payload[num_end..num_end + 2].into();

                        // Optionally, parse the number.
                        let length: usize = payload[num_start..num_end]
                            .parse()
                            .unwrap_or_else(|_| {
                                panic!(
                                    "Failed to convert the splicing size number to usize for the read ID {} given the cs tag: {}.",
                                    read_id, cs_tag
                                )
                            });

                        if reference_strand == Strand::Reverse {
                            donor_splice_site_signal = reverse_complement(&*donor_splice_site_signal);
                            acceptor_splice_site_signal = reverse_complement(&*acceptor_splice_site_signal);
                        }
                        
                        if let Some(left_flank) = last_aligned_read_position {
                            pending_events.push((AlignmentModelEventKind::Splicing, left_flank));
                        }

                        reference_position += length as isize;
                    }
                }
            }

            // Anything still pending is a gap the alignment never closed.
            // It would have to end on a deletion or an intron.
            pending_events.clear();
        }
    }

    /// Close every held event on `right_flank` and file it.
    fn flush_pending_events(
        events: &mut HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,
        events_index: &mut HashMap<ReadPosition, ReadPosition>,
        pending_events: &mut Vec<(AlignmentModelEventKind, ReadPosition)>,
        reference_strand: &Strand,
        right_flank: ReadPosition
    ) {
        if pending_events.is_empty() {
            return;
        }

        for (kind, left_flank) in pending_events.drain(..) {
            let alignment_event: AlignmentModelEvent = if *reference_strand == Strand::Forward {
                AlignmentModelEvent::new(
                    kind,
                    left_flank,
                    right_flank,
                    GraphOperationType::Downstream,
                    GraphOperationType::Upstream
                )
            } else {
                AlignmentModelEvent::new(
                    kind,
                    right_flank,
                    left_flank,
                    GraphOperationType::Upstream,
                    GraphOperationType::Downstream
                )
            };

            Self::add_event(events, events_index, alignment_event);
        }
    }

    /// Add breakpoint events for terminal soft-clips and inter-alignment breakpoints.
    fn apply_breakpoints(
        read_id: ReadID,
        bases: &mut Vec<AlignmentModelBase>,
        events: &mut HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,
        events_index: &mut HashMap<ReadPosition, ReadPosition>,
        alignment_records: &[AlignmentRecord],
        read_sequence: &str
    ) {
        let last_index: usize = alignment_records
            .iter()
            .enumerate()
            .max_by_key(|(_index, alignment_record)| alignment_record.read_end)
            .map(|(index, _alignment_record)| index)
            .unwrap();

        let mut prev_alignment_record: &AlignmentRecord = &alignment_records[0];
        for (i, curr_alignment_record) in alignment_records.iter().enumerate() {
            let reference_chromosome_id: ReferenceChromosomeID = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as ReferenceChromosomeID;
            let reference_strand: Strand = get_alignment_strand(&curr_alignment_record.record);
            let mapping_quality: MappingQuality = get_alignment_mapping_quality(&curr_alignment_record.record);

            // Identify soft-clipped insertion in the first alignment.
            if (i == 0) && (curr_alignment_record.read_start != 0) {
                let expected_softclip_length: u32 = curr_alignment_record.read_start;
                let (softclip_length, reference_position): (u32, ReferencePosition) = if reference_strand == Strand::Reverse {
                    let (is_softclipped, softclip_length): (bool, u32) = get_clipping(&curr_alignment_record.record, false);
                    if !is_softclipped {
                        panic!("The 3' end of the first alignment of read ID {read_id} is expected to be soft-clipped.");
                    }
                    assert_eq!(
                        softclip_length, expected_softclip_length, 
                        "Read start position is expected to be the same as the number of soft-clipped bases."
                    );
                    let reference_position: ReferencePosition = get_alignment_end_position(&curr_alignment_record.record);
                    (softclip_length, reference_position)
                } else {
                    let (is_softclipped, softclip_length): (bool, u32) = get_clipping(&curr_alignment_record.record, true);
                    if !is_softclipped {
                        panic!("The 5' end of the first alignment of read ID {read_id} is expected to be soft-clipped.");
                    }
                    assert_eq!(
                        softclip_length, expected_softclip_length, 
                        "Read start position is expected to be the same as the number of soft-clipped bases."
                    );
                    let reference_position: ReferencePosition = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    (softclip_length, reference_position)
                };
                for j in 0..softclip_length {
                    Self::place_base(
                        bases,
                        j,
                        AlignmentModelBaseKind::Softclip,
                        reference_chromosome_id,
                        reference_position,
                        reference_strand.clone(),
                        mapping_quality,
                        None
                    );
                }
            }

            // Identify soft-clipped insertion in the last alignment.
            if (i == last_index) && ((curr_alignment_record.read_end) != read_sequence.len() as u32 - 1) {
                let expected_softclip_length: u32 = read_sequence.len() as u32 - curr_alignment_record.read_end - 1;
                let (softclip_length, reference_position): (u32, ReferencePosition) = if reference_strand == Strand::Reverse {
                    let (is_softclipped, softclip_length): (bool, u32) = get_clipping(&curr_alignment_record.record, true);
                    if !is_softclipped {
                        panic!("The 5' end of the last alignment of read ID {read_id} is expected to be soft-clipped.");
                    }
                    assert_eq!(
                        softclip_length, expected_softclip_length, 
                        "(Read length - alignment's last read position - 1) is expected to match the number of soft-clipped bases."
                    );
                    let reference_position: ReferencePosition = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    (softclip_length, reference_position)
                } else {
                    let (is_softclipped, softclip_length): (bool, u32) = get_clipping(&curr_alignment_record.record, false);
                    if !is_softclipped {
                        panic!("The 3' end of the last alignment of read ID {read_id} is expected to be soft-clipped.");
                    }
                    assert_eq!(
                        softclip_length, expected_softclip_length, 
                        "(Read length - alignment's last read position - 1) is expected to match the number of soft-clipped bases."
                    );
                    let reference_position: ReferencePosition = get_alignment_end_position(&curr_alignment_record.record);
                    (softclip_length, reference_position)
                };
                let start: ReadPosition = read_sequence.len() as u32 - softclip_length;
                let end: ReadPosition = read_sequence.len() as u32;
                for j in start..end {
                    Self::place_base(
                        bases,
                        j,
                        AlignmentModelBaseKind::Softclip,
                        reference_chromosome_id,
                        reference_position,
                        reference_strand.clone(),
                        mapping_quality,
                        None
                    );
                }
            }

            // Identify breakpoints (soft-clipping) between alignments.
            if i > 0 {
                let mut bnd_1_read_position: ReadPosition = prev_alignment_record.read_end;
                let mut bnd_2_read_position: ReadPosition = curr_alignment_record.read_start;

                // Check if curr_alignment_record is completely contained within prev_alignment_record.
                let contained_1: bool = interval_contains(
                    prev_alignment_record.read_start,
                    prev_alignment_record.read_end,
                    curr_alignment_record.read_start,
                    curr_alignment_record.read_end
                );
                let contained_2: bool = interval_contains(
                    curr_alignment_record.read_start,
                    curr_alignment_record.read_end,
                    prev_alignment_record.read_start,
                    prev_alignment_record.read_end
                );

                if contained_1 == false && contained_2 == false {
                    // Check if the previous and the current alignments overlap.
                    let alignments_overlap: bool = overlaps(
                        prev_alignment_record.read_start as isize,
                        prev_alignment_record.read_end as isize,
                        curr_alignment_record.read_start as isize,
                        curr_alignment_record.read_end as isize
                    );
                    let mut insertion: Box<str> = "".to_string().into_boxed_str();
                    if alignments_overlap {
                        // If the previous and the current alignment records overlap,
                        // treat the overlapping part as an insertion.
                        let (overlap_start, overlap_end) = find_overlap(
                            (prev_alignment_record.read_start as isize, prev_alignment_record.read_end as isize),
                            (curr_alignment_record.read_start as isize, curr_alignment_record.read_end as isize)
                        ).unwrap();
                        insertion = read_sequence[(overlap_start as usize)..=(overlap_end as usize)].to_string().into_boxed_str();

                        // Retreat read positions by the length of the insertion and mark each an insertion.
                        for j in overlap_start..=overlap_end {
                            Self::set_base_kind(bases, j as ReadPosition, AlignmentModelBaseKind::Softclip);
                            bnd_1_read_position -= 1;
                            bnd_2_read_position += 1;
                        }
                    } else {
                        // Check if an insertion (i.e. unaligned bases) exists between the breakpoints.
                        if prev_alignment_record.read_end + 1 != curr_alignment_record.read_start &&
                            prev_alignment_record.read_end < curr_alignment_record.read_start {
                            insertion = read_sequence[
                                (prev_alignment_record.read_end as usize) + 1..=(curr_alignment_record.read_start as usize) - 1
                            ].to_string().into_boxed_str();

                            // Mark each unaligned base an insertion.
                            for j in (prev_alignment_record.read_end) + 1..=(curr_alignment_record.read_start) - 1 {
                                Self::set_base_kind(bases, j, AlignmentModelBaseKind::Softclip);
                            }
                        }
                    }

                    assert!(
                        bnd_1_read_position < bnd_2_read_position, 
                        "{} < {} not satisfied", bnd_1_read_position, bnd_2_read_position
                    );
                    assert!(
                        bnd_2_read_position < read_sequence.len() as u32, 
                        "{} < {} not satisfied", bnd_2_read_position, read_sequence.len()
                    );

                    let bnd_1_strand: Strand = match bases[bnd_1_read_position as usize].get_placement() {
                        AlignmentModelBasePlacement::Placed { strand, .. } => strand.clone(),
                        AlignmentModelBasePlacement::Unplaced => panic!(
                            "read {}: breakpoint base at read position {} is not placed on the reference",
                            read_id, bnd_1_read_position
                        )
                    };
                    let bnd_2_strand: Strand = match bases[bnd_2_read_position as usize].get_placement() {
                        AlignmentModelBasePlacement::Placed { strand, .. } => strand.clone(),
                        AlignmentModelBasePlacement::Unplaced => panic!(
                            "read {}: breakpoint base at read position {} is not placed on the reference",
                            read_id, bnd_2_read_position
                        )
                    };

                    let bnd_1_operation: GraphOperationType = GraphOperationType::for_breakpoint(
                        &bnd_1_strand,
                        false
                    );
                    let bnd_2_operation: GraphOperationType = GraphOperationType::for_breakpoint(
                        &bnd_2_strand,
                        true
                    );

                    let alignment_event: AlignmentModelEvent = AlignmentModelEvent::new(
                        AlignmentModelEventKind::Breakpoint,
                        bnd_1_read_position,
                        bnd_2_read_position,
                        bnd_1_operation.clone(),
                        bnd_2_operation.clone()
                    );

                    Self::add_event(events, events_index, alignment_event);
                }
            }

            // A record inside the read span of the one before adds no read bases: the next record
            // joins the one before it.
            let is_contained: bool = interval_contains(
                prev_alignment_record.read_start,
                prev_alignment_record.read_end,
                curr_alignment_record.read_start,
                curr_alignment_record.read_end
            );
            if !is_contained {
                prev_alignment_record = curr_alignment_record;
            }
        }
    }

    /// Set a base's kind and place it on the reference.
    fn place_base(
        bases: &mut Vec<AlignmentModelBase>,
        read_position: ReadPosition,
        kind: AlignmentModelBaseKind,
        chromosome_id: ReferenceChromosomeID,
        position: ReferencePosition,
        strand: Strand,
        mapping_quality: MappingQuality,
        reference_nucleotide: Option<Nucleotide>
    ) {
        let base: &mut AlignmentModelBase = &mut bases[read_position as usize];
        base.set_kind(kind);
        base.place(chromosome_id, position, strand, mapping_quality);
        base.set_reference_nucleotide(reference_nucleotide);
        base.set_deleted_reference_bases(Vec::new());
    }

    /// Set a base's kind without changing its placement.
    fn set_base_kind(
        bases: &mut Vec<AlignmentModelBase>,
        read_position: ReadPosition,
        kind: AlignmentModelBaseKind
    ) {
        bases[read_position as usize].set_kind(kind);
    }

    /// Drop every event whose flanking read bases are about to be re-placed.
    fn purge_events(
        events: &mut HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,
        events_index: &mut HashMap<ReadPosition, ReadPosition>,
        read_start: ReadPosition,
        read_end: ReadPosition
    ) {
        if events.is_empty() {
            return;
        }

        let stale: Vec<(ReadPosition, ReadPosition)> = events
            .keys()
            .filter(|(read_position_1, read_position_2)| {
                (*read_position_1 >= read_start && *read_position_1 <= read_end) ||
                    (*read_position_2 >= read_start && *read_position_2 <= read_end)
            })
            .copied()
            .collect();

        if stale.is_empty() {
            return;
        }

        for key in stale.iter() {
            events.remove(key);

            if events_index.get(&key.0) == Some(&key.1) {
                events_index.remove(&key.0);
            }
            if events_index.get(&key.1) == Some(&key.0) {
                events_index.remove(&key.1);
            }
        }

        // Refill the slots the purge emptied, leaving live entries alone. Sorted so that two
        // survivors competing for one freed slot resolve the same way on every run.
        let mut surviving: Vec<(ReadPosition, ReadPosition)> = events.keys().copied().collect();
        surviving.sort_unstable();
        for (read_position_1, read_position_2) in surviving.iter() {
            events_index.entry(*read_position_1).or_insert(*read_position_2);
            events_index.entry(*read_position_2).or_insert(*read_position_1);
        }
    }

    /// Insert an event, maintaining the bidirectional `events_index` and the
    /// splicing-supersedes-deletion precedence rule.
    fn add_event(
        events: &mut HashMap<(ReadPosition, ReadPosition), AlignmentModelEvent>,
        events_index: &mut HashMap<ReadPosition, ReadPosition>,
        event: AlignmentModelEvent
    ) {
        let key: (ReadPosition, ReadPosition) = (event.get_prev_read_position(), event.get_next_read_position());

        if let Some(existing) = events.get(&key) {
            let incoming_is_splice: bool =
                *event.get_kind() == AlignmentModelEventKind::Splicing;
            let existing_is_splice: bool =
                *existing.get_kind() == AlignmentModelEventKind::Splicing;

            // Splicing supersedes deletion; deletion never supersedes splicing.
            if existing_is_splice && !incoming_is_splice {
                return;
            }
        }

        events_index.insert(key.0, key.1);
        events_index.insert(key.1, key.0);
        events.insert(key, event);
    }
}


fn get_clipping(record: &bam::Record, left: bool) -> (bool, u32) {
    let operations: Vec<(Kind, u32)> = record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len() as u32)).unwrap())
        .collect();
    let is_clip = |&&(kind, _): &&(Kind, u32)| matches!(kind, Kind::SoftClip | Kind::HardClip);
    let length: u32 = if left {
        operations.iter().take_while(is_clip).map(|(_, len)| len).sum()
    } else {
        operations.iter().rev().take_while(is_clip).map(|(_, len)| len).sum()
    };
    (length > 0, length)
}
