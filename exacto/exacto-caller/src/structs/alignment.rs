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


use bstr::ByteSlice;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_sam::alignment::Record;
use regex::Regex;
use std::str::FromStr;

use crate::prelude::*;


/// Represents the alignment of a single read against a reference.
///
/// This struct stores both the original FASTX information (sequence and
/// quality scores) as well as processed alignment results. It is designed
/// to keep traceability from the raw input read to its aligned form.
///
/// # Fields
/// * `read_id` - A unique identifier for the read.
///
/// * `read_sequence` - The original nucleotide sequence from the FASTX read.
///   Stored as a `Box<str>` to reduce heap allocation overhead.
///
/// * `base_quality_scores` - The per-base quality scores associated with the
///   original read sequence.
///
/// * `alignment_records` - A collection of `AlignmentRecord` objects, each
///   representing one alignment placement of the read (e.g., primary or
///   secondary/supplementary alignments).
///
/// * `alignment_structure` - The consensus aligned structure representation,
///   which may include gaps, mismatches, or other transformations introduced
///   during alignment.
#[derive(Debug)]
pub struct Alignment {
    read_id: usize,
    read_sequence: Box<str>,
    base_quality_scores: Vec<u8>,
    alignment_records: Vec<AlignmentRecord>,
    alignment_structure: AlignmentStructure
}

/// API methods
impl Alignment {

    pub fn new(
        read_id: usize,
        read_sequence: &str,
        base_quality_scores: &Vec<u8>,
        records: &Vec<bam::Record>
    ) -> Self {
        assert!(!records.is_empty());

        // Step 1. Make sure all the BAM records come from the same read ID
        let first_read_id = std::str::from_utf8(records[0].name().unwrap().as_bytes()).unwrap();
        for record in records.iter().skip(1) {
            let read_id = std::str::from_utf8(record.name().unwrap().as_bytes()).unwrap();
            assert_eq!(read_id, first_read_id, "Not all records have the same read ID.");
        }

        // Step 2. Build alignment records
        let alignment_records: Vec<AlignmentRecord> = Self::build_alignment_records(
            read_sequence,
            records
        );

        // Step 3. Initialize alignment sequence
        let mut alignment_structure: AlignmentStructure = Self::init_alignment_structure(
            read_id, 
            read_sequence,
            base_quality_scores
        );
        
        // Step 4. Identify local variants (using CS tags)
        Self::identify_local_variants(
            &mut alignment_structure, 
            &alignment_records
        );

        // Step 3. Identify breakpoint variants (using SA tags and soft-clipped bases)
        Self::identify_breakpoint_variants(
            &mut alignment_structure,
            &alignment_records,
            read_id,
            read_sequence
        );

        // Step 4. Check if the read sequence matches the original read sequence
        assert_eq!(
            read_sequence.to_uppercase(),
            alignment_structure.get_read_sequence().to_uppercase(),
            "read_sequence:\n{}\nalignment_structure.get_read_sequence():\n{}",
            read_sequence.to_uppercase(),
            alignment_structure.get_read_sequence().to_uppercase()
        );
        
        Self {
            read_id,
            read_sequence: read_sequence.into(),
            base_quality_scores: base_quality_scores.clone(),
            alignment_records,
            alignment_structure: alignment_structure
        }
    }
    
    pub fn get_alignment_records(&self) -> &Vec<AlignmentRecord> {
        &self.alignment_records
    }
    
    pub fn get_alignment_records_count(&self) -> usize {
        self.alignment_records.len()
    }

    pub fn get_alignment_structure(&self) -> &AlignmentStructure {
        &self.alignment_structure
    }
    
    pub fn get_base_quality_scores(&self) -> &Vec<u8> {
        &self.base_quality_scores
    }

    pub fn get_read_id(&self) -> usize {
        self.read_id
    }
    
    pub fn get_read_length(&self) -> usize {
        self.read_sequence.len()
    }
    
    pub fn get_read_sequence(&self) -> &str {
        &*self.read_sequence
    }
}

// Helper methods
impl Alignment {
    fn build_alignment_records(
        read_sequence: &str,
        records: &Vec<bam::Record>
    ) -> Vec<AlignmentRecord> {
        let mut alignment_records: Vec<AlignmentRecord> = Vec::new();

        // Step 1. Identify alignment records
        for record in records {
            // Get the aligned sequence
            let aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(&record).into();

            // Identify all start positions between aligned sequence and original read sequence
            let start_positions: Vec<usize> = find_substring_positions(&*read_sequence.to_uppercase(), &*aligned_sequence.to_uppercase());
            assert!(!start_positions.is_empty(), "Could not find the aligned sequence in the original read sequence.");

            // Get left and right soft-clipping information of the current record
            let left_softclipping: (bool, usize) = get_left_softclipping(&record);
            let right_softclipping: (bool, usize) = get_right_softclipping(&record);

            // If there are multiple start positions, find where the aligned sequence starts
            // on the original read sequence
            let mut read_start: usize = 0;
            let mut read_end: usize = 0;
            let reference_strand: Strand = get_alignment_strand(&record);
            for start_position in start_positions.iter() {
                // Check if the current start position aligns with the current alignment record
                let end_position: usize = *start_position + aligned_sequence.len() - 1;
                let num_left_bases: usize = *start_position;
                let num_right_bases: usize = read_sequence.len() - end_position - 1;

                let mut aligned: bool = true;
                if reference_strand == Strand::Reverse {
                    if (left_softclipping.0 && left_softclipping.1 != num_right_bases) ||
                        (right_softclipping.0 && right_softclipping.1 != num_left_bases) {
                        aligned = false;
                    }
                } else {
                    if (left_softclipping.0 && left_softclipping.1 != num_left_bases) ||
                        (right_softclipping.0 && right_softclipping.1 != num_right_bases) {
                        aligned = false;
                    }
                }

                if aligned {
                    read_start = *start_position as usize;
                    read_end = read_start + (aligned_sequence.len() as usize) - 1;
                    break;
                }
            }

            assert!(read_start != read_end, "read_start should not be the same as read_end.");
            assert!(&*aligned_sequence == read_sequence[(read_start as usize)..(read_end as usize)+1].to_string(), "Aligned sequence does not match the identified part of the original read sequence.");

            let alignment_record: AlignmentRecord = AlignmentRecord::new(
                read_start,
                read_end,
                reference_strand,
                record.clone()
            );

            alignment_records.push(alignment_record);
        }

        // Step 2. Sort alignment records by read start position
        alignment_records.sort_by_key(|alignment| alignment.read_start);

        alignment_records
    }

    fn init_alignment_structure(
        read_id: usize,
        read_sequence: &str,
        base_quality_scores: &Vec<u8>
    ) -> AlignmentStructure {
        let mut alignment_structure: AlignmentStructure = AlignmentStructure::new(read_id);
        for (i, s) in read_sequence.chars().enumerate() {
            let nucleotide: Nucleotide = Nucleotide::from_str(s.to_string().as_str()).unwrap();
            let base_quality: u8 = base_quality_scores[i];
            let alignment_base: AlignmentStructureBase = AlignmentStructureBase::new(
                i as usize,
                nucleotide,
                base_quality
            );
            alignment_structure.add_base(alignment_base);
        }
        alignment_structure
    }

    /// Identify breakpoint variants using softclipped bases.
    fn identify_breakpoint_variants(
        alignment_structure: &mut AlignmentStructure,
        alignment_records: &Vec<AlignmentRecord>,
        read_id: usize,
        read_sequence: &str
    ) {
        let mut prev_alignment_record: &AlignmentRecord = &alignment_records[0];
        for (i, curr_alignment_record) in alignment_records.iter().enumerate() {
            let reference_chromosome_id: u16 = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
            let reference_strand: Strand = get_alignment_strand(&curr_alignment_record.record);
            let mapping_quality: usize = get_alignment_mapping_quality(&curr_alignment_record.record);

            // Identify soft-clipped insertion in the first alignment
            if (i == 0) && (curr_alignment_record.read_start != 0) {
                let expected_softclip_length: usize = curr_alignment_record.read_start as usize;
                let (is_softclipped, softclip_length, reference_position) = if reference_strand == Strand::Reverse {
                    let (is_softclipped, softclip_length) = get_right_softclipping(&curr_alignment_record.record);
                    assert!(is_softclipped, "The 3' end of the first alignment (read ID: {}) is soft-clipped.", read_id);
                    assert_eq!(softclip_length, expected_softclip_length, "Read start position is expected to be the same as the number of soft-clipped bases.");
                    let reference_position: usize = get_alignment_end_position(&curr_alignment_record.record);
                    (is_softclipped, softclip_length, reference_position)
                } else {
                    let (is_softclipped, softclip_length) = get_left_softclipping(&curr_alignment_record.record);
                    assert!(is_softclipped, "The 5' end of the first alignment (read ID: {}) is soft-clipped.", read_id);
                    assert_eq!(softclip_length, expected_softclip_length, "Read start position is expected to be the same as the number of soft-clipped bases.");
                    let reference_position: usize = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    (is_softclipped, softclip_length, reference_position)
                };
                for j in 0..softclip_length {
                    let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(j as usize);
                    alignment_base.set_mapping_quality(mapping_quality);
                    alignment_base.set_reference_chromosome_id(reference_chromosome_id);
                    alignment_base.set_reference_position(reference_position);
                    alignment_base.set_reference_strand(reference_strand.clone());
                    alignment_base.set_is_soft_clipped(true);
                    alignment_base.set_kind(AlignmentStructureBaseKind::Insertion);
                }
            }

            // Identify soft-clipped insertion in the last alignment
            if (i == alignment_records.len() - 1) && ((curr_alignment_record.read_end as usize) != read_sequence.len() - 1) {
                let expected_softclip_length: usize = read_sequence.len() - (curr_alignment_record.read_end as usize) - 1;
                let (is_softclipped, softclip_length, reference_position) = if reference_strand == Strand::Reverse {
                    let (is_softclipped, softclip_length) = get_left_softclipping(&curr_alignment_record.record);
                    assert!(is_softclipped, "The 5' end of the last alignment (read ID: {}) is soft-clipped.", read_id);
                    assert_eq!(softclip_length, expected_softclip_length, "(Read length - alignment's last read position - 1) is expected to match the number of soft-clipped bases.");
                    let reference_position: usize = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    (is_softclipped, softclip_length, reference_position)
                } else {
                    let (is_softclipped, softclip_length) = get_right_softclipping(&curr_alignment_record.record);
                    assert!(is_softclipped, "The 3' end of the last alignment (read ID: {}) is soft-clipped.", read_id);
                    assert_eq!(softclip_length, expected_softclip_length, "(Read length - alignment's last read position - 1) is expected to match the number of soft-clipped bases.");
                    let reference_position = get_alignment_end_position(&curr_alignment_record.record);
                    (is_softclipped, softclip_length, reference_position)
                };
                let start: usize = read_sequence.len() - softclip_length;
                let end: usize = read_sequence.len();
                for j in start..end {
                    let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(j as usize);
                    alignment_base.set_mapping_quality(mapping_quality);
                    alignment_base.set_reference_chromosome_id(reference_chromosome_id);
                    alignment_base.set_reference_position(reference_position);
                    alignment_base.set_reference_strand(reference_strand.clone());
                    alignment_base.set_is_soft_clipped(true);
                    alignment_base.set_kind(AlignmentStructureBaseKind::Insertion);
                }
            }

            // Identify breakpoints (soft-clipping) between alignments
            if i > 0 {
                let mut bnd_1_read_position: usize = prev_alignment_record.read_end;
                let mut bnd_2_read_position: usize = curr_alignment_record.read_start;

                // Check if the previous and the current alignments overlap
                let alignments_overlap: bool = overlaps(
                    prev_alignment_record.read_start as isize,
                    prev_alignment_record.read_end as isize,
                    curr_alignment_record.read_start as isize,
                    curr_alignment_record.read_end as isize
                );
                if alignments_overlap {
                    // If the previous and the current alignment records overlap,
                    // treat the overlapping part as an insertion
                    let (overlap_start,overlap_end) = find_overlap(
                        (prev_alignment_record.read_start as isize,prev_alignment_record.read_end as isize),
                        (curr_alignment_record.read_start as isize,curr_alignment_record.read_end as isize)
                    ).unwrap();
                    let insertion: Box<str> = read_sequence[(overlap_start as usize)..=(overlap_end as usize)].to_string().into_boxed_str();

                    // Retreat read positions by the length of the insertion and mark each an insertion
                    for j in overlap_start..=overlap_end {
                        let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(j as usize);
                        alignment_base.set_is_embedded_insertion(true);
                        alignment_base.set_kind(AlignmentStructureBaseKind::Insertion);
                        bnd_1_read_position -= 1;
                        bnd_2_read_position += 1;
                    }
                } else {
                    // Check if an insertion (i.e. unaligned bases) exists between the breakpoints
                    if prev_alignment_record.read_end + 1 != curr_alignment_record.read_start &&
                        prev_alignment_record.read_end < curr_alignment_record.read_start {
                        let insertion: Box<str> = read_sequence[(prev_alignment_record.read_end as usize) + 1..=(curr_alignment_record.read_start as usize) - 1].to_string().into_boxed_str();

                        // Mark each unaligned base an insertion
                        for j in (prev_alignment_record.read_end as usize) + 1..=(curr_alignment_record.read_start as usize) - 1 {
                            let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(j as usize);
                            alignment_base.set_is_embedded_insertion(true);
                            alignment_base.set_kind(AlignmentStructureBaseKind::Insertion);
                        }
                    }
                }

                assert!(bnd_1_read_position < bnd_2_read_position);

                let bnd_1_base: &AlignmentStructureBase = alignment_structure.get_base(bnd_1_read_position);
                let bnd_2_base: &AlignmentStructureBase = alignment_structure.get_base(bnd_2_read_position);
                let bnd_1_operation: GraphOperationType = match bnd_1_base.get_reference_strand().as_ref().unwrap() {
                    Strand::Forward => GraphOperationType::Downstream,
                    Strand::Reverse => GraphOperationType::Upstream,
                    Strand::Both => panic!("Unexpected strand: {}", bnd_1_base.get_reference_strand().as_ref().unwrap().as_str()),
                    Strand::Unknown => panic!("Unexpected strand: {}", bnd_1_base.get_reference_strand().as_ref().unwrap().as_str())
                };
                let bnd_2_operation: GraphOperationType = match bnd_2_base.get_reference_strand().as_ref().unwrap() {
                    Strand::Forward => GraphOperationType::Upstream,
                    Strand::Reverse => GraphOperationType::Downstream,
                    Strand::Both => panic!("Unexpected strand: {}", bnd_2_base.get_reference_strand().as_ref().unwrap().as_str()),
                    Strand::Unknown => panic!("Unexpected strand: {}", bnd_2_base.get_reference_strand().as_ref().unwrap().as_str())
                };

                let alignment_event: AlignmentStructureEvent = AlignmentStructureEvent::new(
                    AlignmentStructureEventKind::Breakpoint,
                    bnd_1_read_position,
                    bnd_2_read_position,
                    bnd_1_operation.clone(),
                    bnd_2_operation.clone()
                );

                alignment_structure.add_event(alignment_event);
            }
            prev_alignment_record = curr_alignment_record;
        }
    }
    
    /// Identify local variants (SNVs, insertions, splicing, and deletions).
    fn identify_local_variants(
        alignment_structure: &mut AlignmentStructure,
        alignment_records: &Vec<AlignmentRecord>
    ) {
        // Step 1. Update alignment sequence
        for alignment_record in alignment_records.iter() {
            let mut reference_position: isize = get_alignment_start_position(&alignment_record.record) as isize - 1;
            let reference_chromosome_id: u16 = alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
            let reference_strand: Strand = get_alignment_strand(&alignment_record.record);
            let mapping_quality: usize = get_alignment_mapping_quality(&alignment_record.record);
            let cs_tag: String = get_tag_value(&alignment_record.record, "cs")
                .expect("Could not find the CS tag.")
                .to_string();
            let mut read_position: isize = if reference_strand == Strand::Forward {
                alignment_record.read_start as isize - 1
            } else {
                alignment_record.read_end as isize + 1
            };

            // Identify SNVs, insertions, deletions, and splicing in the CS tag
            let re = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
            for cap in re.captures_iter(&cs_tag) {
                let token = &cap[0];
                let mut chars = token.chars();
                let cs_tag_kind: CSTagKind = CSTagKind::from_str(chars.next().unwrap().to_string().as_str()).unwrap();
                let payload = chars.as_str();

                match cs_tag_kind {
                    CSTagKind::Match => {
                        let length: isize = payload.parse().unwrap();
                        for _ in 0..length {
                            read_position += if reference_strand == Strand::Forward { 1 } else { -1 };
                            reference_position += 1;
                            let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(read_position as usize);
                            alignment_base.set_mapping_quality(mapping_quality);
                            alignment_base.set_reference_chromosome_id(reference_chromosome_id);
                            alignment_base.set_reference_position(reference_position as usize);
                            alignment_base.set_reference_strand(reference_strand.clone());
                            alignment_base.set_kind(AlignmentStructureBaseKind::Match);
                        }
                    },
                    CSTagKind::Mismatch => {
                        let alleles: Vec<char> = payload.chars().collect();
                        assert_eq!(alleles.len(), 2, "1 reference allele and 1 alternate allele expected.");
                        read_position += if reference_strand == Strand::Forward { 1 } else { -1 };
                        reference_position += 1;
                        let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(read_position as usize);
                        alignment_base.set_mapping_quality(mapping_quality);
                        alignment_base.set_reference_chromosome_id(reference_chromosome_id);
                        alignment_base.set_reference_position(reference_position as usize);
                        alignment_base.set_reference_strand(reference_strand.clone());
                        alignment_base.set_kind(AlignmentStructureBaseKind::Mismatch);
                    },
                    CSTagKind::Insertion => {
                        let insertion: String = payload.to_string();
                        let length: usize = insertion.chars().count();
                        for _ in 0..length {
                            read_position += if reference_strand == Strand::Forward { 1 } else { -1 };
                            let alignment_base: &mut AlignmentStructureBase = alignment_structure.get_mut_base(read_position as usize);
                            alignment_base.set_mapping_quality(mapping_quality);
                            alignment_base.set_reference_chromosome_id(reference_chromosome_id);
                            alignment_base.set_reference_position(reference_position as usize);
                            alignment_base.set_reference_strand(reference_strand.clone());
                            alignment_base.set_kind(AlignmentStructureBaseKind::Insertion);
                        }
                    }
                    CSTagKind::Deletion => {
                        let length: usize = payload.chars().count();
                        let read_position_1: usize = read_position as usize;
                        let read_position_2: usize = if reference_strand == Strand::Forward { read_position as usize + 1 } else { read_position as usize - 1 };
                        
                        let alignment_event: AlignmentStructureEvent = if reference_strand == Strand::Forward {
                            AlignmentStructureEvent::new(
                                AlignmentStructureEventKind::Deletion,
                                read_position_1,
                                read_position_2,
                                GraphOperationType::Downstream,
                                GraphOperationType::Upstream
                            )
                        } else {
                            AlignmentStructureEvent::new(
                                AlignmentStructureEventKind::Deletion,
                                read_position_2,
                                read_position_1,
                                GraphOperationType::Upstream,
                                GraphOperationType::Downstream
                            )
                        };

                        alignment_structure.add_event(alignment_event);

                        reference_position += length as isize;
                    }
                    CSTagKind::Splicing => {
                        let re_splicing = Regex::new(r"\d+").unwrap();
                        let caps = re_splicing.find(&payload).expect("No numerical value found");

                        let num_start = caps.start();
                        let num_end = caps.end();

                        // Extract donor splice site signal (2 letters before the number)
                        let mut donor_splice_site_signal: Box<str> = payload[num_start - 2..num_start].into();

                        // Extract acceptor splice site signal (2 letters after the number)
                        let mut acceptor_splice_site_signal: Box<str> = payload[num_end..num_end + 2].into();

                        // Optionally, parse the number
                        let length: usize = payload[num_start..num_end]
                            .parse()
                            .expect("Failed to convert the splicing size number to usize");

                        if reference_strand == Strand::Reverse {
                            donor_splice_site_signal = reverse_complement(&*donor_splice_site_signal);
                            acceptor_splice_site_signal = reverse_complement(&*acceptor_splice_site_signal);
                        }

                        let read_position_1: usize = read_position as usize;
                        let read_position_2: usize = if reference_strand == Strand::Forward { read_position as usize + 1 } else { read_position as usize - 1 };

                        let alignment_event: AlignmentStructureEvent = if reference_strand == Strand::Forward {
                            AlignmentStructureEvent::new(
                                AlignmentStructureEventKind::Splicing,
                                read_position_1,
                                read_position_2,
                                GraphOperationType::Downstream,
                                GraphOperationType::Upstream
                            )
                        } else {
                            AlignmentStructureEvent::new(
                                AlignmentStructureEventKind::Splicing,
                                read_position_2,
                                read_position_1,
                                GraphOperationType::Upstream,
                                GraphOperationType::Downstream
                            )
                        };

                        alignment_structure.add_event(alignment_event);

                        reference_position += length as isize;
                    }
                }
            }
        }
    }
}
