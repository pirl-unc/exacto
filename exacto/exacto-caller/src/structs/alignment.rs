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
use bstr::ByteSlice;
use exacto_util::prelude::*;
use noodles_bam as bam;
use noodles_core::{Region, Position};
use noodles_fasta::io::indexed_reader::Builder;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::record::Flags;
use regex::Regex;
use std::collections::{HashSet, VecDeque};

use crate::algorithms::variant_calling_rna::*;
use crate::common::bam::*;
use crate::common::constants::*;
use crate::prelude::ReferenceTranscriptMatch;
use crate::structs::alignment_record::AlignmentRecord;
use crate::structs::transcript_model_exon::TranscriptModelExon;
use crate::structs::sequence_operation::SequenceOperation;
use crate::structs::transcript_model_splice_junction::TranscriptModelSpliceJunction;
use crate::structs::variant_record::VariantRecord;


#[derive(Debug)]
pub struct Alignment {
    pub read_id: usize,
    pub original_read_sequence: Box<str>,
    pub quality_scores: Vec<u8>,
    pub alignment_records: Vec<AlignmentRecord>
}

impl Alignment {
    pub fn new(
        read_id: usize,
        original_read_sequence: Box<str>,
        quality_scores: Vec<u8>,
        records: Vec<bam::Record>
    ) -> Self {
        assert!(!records.is_empty());

        // Step 1. Make sure all the BAM records come from the same read ID
        let first_read_id: &str = std::str::from_utf8(records[0].name().unwrap().as_bytes()).unwrap();
        for record in records.iter().skip(1) {
            assert!(std::str::from_utf8(record.name().unwrap().as_bytes()).unwrap() == first_read_id, "Not all records have the same read ID.");
        }

        // Step 2. Identify alignment records
        let mut alignment_records: Vec<AlignmentRecord> = Vec::new();
        for record in records {
            let mut aligned_sequence: Box<str> = get_aligned_sequence_from_cigar(&record).to_uppercase().into_boxed_str();
            let mut read_start: usize = 0;
            let mut read_end: usize = 0;
            let mut reverse_complemented: bool = false;
            let start_positions: Vec<usize> = find_substring_positions(&*original_read_sequence, &*aligned_sequence);

            assert!(!start_positions.is_empty(), "Could not find the aligned sequence in the original read sequence.");

            for start_position in start_positions.iter() {
                // Check if the current start position aligns with the current alignment record
                let end_position: usize = *start_position + aligned_sequence.len() - 1;
                let num_left_bases: usize = *start_position;
                let num_right_bases: usize = original_read_sequence.len() - end_position - 1;
                let left_softclipping: (bool,usize) = get_left_softclipping(&record);
                let right_softclipping: (bool,usize) = get_right_softclipping(&record);

                let mut aligned: bool = true;
                if record.flags().is_reverse_complemented() {
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
                    read_start = *start_position;
                    read_end = read_start + aligned_sequence.len() - 1;
                    reverse_complemented = record.flags().is_reverse_complemented();
                    break;
                }
            }

            assert!(read_start != read_end, "read_start should not be the same as read_end.");

            let alignment_record: AlignmentRecord = AlignmentRecord::new(
                read_start,
                read_end,
                reverse_complemented,
                record
            );

            alignment_records.push(alignment_record);
        }

        // Step 3. Sort the alignment records
        alignment_records.sort_by_key(|alignment| alignment.read_start);

        Self {
            read_id,
            original_read_sequence,
            quality_scores,
            alignment_records: alignment_records
        }
    }

    pub fn get_alignment_records_count(&self) -> usize {
        self.alignment_records.len()
    }

    pub fn get_read_length(&self) -> usize {
        self.original_read_sequence.len()
    }

    /// Identify exonic boundaries.
    pub fn identify_exons(
        &self,
        min_mapping_quality: usize
    ) -> Vec<TranscriptModelExon> {
        // Step 1. Identify exonic boundaries
        let mut exons_: VecDeque<TranscriptModelExon> = VecDeque::new();
        for alignment_record in self.alignment_records.iter() {
            // Check if the mapping quality meets the minimum mapping quality
            if min_mapping_quality > alignment_record.record.mapping_quality().unwrap().get() as usize {
                continue;
            }

            // Get the chromosome name
            let chromosome_id: u16 = alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;

            // Get the alignment flag
            let mut strand: Strands = Strands::Forward;
            for flag in alignment_record.record.flags() {
                if flag == Flags::REVERSE_COMPLEMENTED {
                    strand = Strands::Reverse;
                    break;
                }
            }

            // Get the CS tag
            let cs_tag: String;
            if let Some(value) = get_tag_value(&alignment_record.record, "cs") {
                cs_tag = value.to_string();
            } else {
                panic!("Could not find the CS tag.");
            }

            // Get the start reference position
            // curr_reference_pos always points at the end of the last variant
            let mut curr_reference_pos: isize = match alignment_record.record.alignment_start().unwrap() {
                Ok(s) => {
                    s.get() as isize
                },
                Err(e) => {
                    panic!("Could not fetch the start position");
                },
            };
            let mut curr_exon_start: usize = curr_reference_pos as usize;
            curr_reference_pos -= 1;

            // Identify exons
            let re = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
            for value in re.captures_iter(&cs_tag) {
                if value[0].contains(":") {
                    // Increment current reference position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut chars = value[0].chars();
                    chars.next();

                    // Increment current position by the number of matched nucleotides
                    let num_matched_nucleotides: isize = chars.as_str().parse::<isize>().unwrap();
                    curr_reference_pos += num_matched_nucleotides - 1;
                } else if value[0].contains("*") {
                    // Increment current reference position
                    curr_reference_pos += 1;
                } else if value[0].contains("+") {
                    // Do nothing
                } else if value[0].contains("-") {
                    // Increment current position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Reference and alternate alleles
                    let sequence: String = value_chars.as_str().to_string().to_uppercase();
                    let variant_size: usize = sequence.chars().count();

                    // Increment current reference position
                    curr_reference_pos += (variant_size - 1) as isize;
                } else if value[0].contains("~") {
                    // Record an exon
                    let curr_exon_end: usize =  curr_reference_pos as usize;
                    let exon: TranscriptModelExon = TranscriptModelExon::new(
                        chromosome_id,
                        curr_exon_start as u32,
                        curr_exon_end as u32,
                        0,
                        strand.clone()
                    );
                    if strand == Strands::Reverse {
                        exons_.push_front(exon);
                    } else {
                        exons_.push_back(exon);
                    }

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Get splicing size
                    let splicing_str: String = value_chars.as_str().to_string();
                    let re_splicing = Regex::new(r"\d+").unwrap(); // r"\d+" matches one or more digits
                    let splicing_size: usize = match re_splicing.find(&splicing_str) {
                        Some(matched) => {
                            let num_str = &splicing_str[matched.start()..matched.end()];
                            let size: usize = match num_str.parse::<usize>() {
                                Ok(num) => num,
                                Err(e) => {
                                    panic!("Failed to convert the splicing size number to usize");
                                },
                            };
                            size
                        },
                        None => {
                            panic!("No numerical value was found in the splicing string: {}", splicing_str);
                        }
                    };

                    // Increment current reference position
                    curr_reference_pos += 1;
                    curr_reference_pos += (splicing_size - 1) as isize;

                    // Update current exon start position
                    curr_exon_start = (curr_reference_pos + 1) as usize;
                } else {
                    panic!("Unknown delimiter for CS tag: {}", value[0].to_string());
                }
            }

            // Include the first/last exon
            let curr_exon_end: usize =  curr_reference_pos as usize;
            let exon: TranscriptModelExon = TranscriptModelExon::new(
                chromosome_id,
                curr_exon_start as u32,
                curr_exon_end as u32,
                0,
                strand.clone()
            );
            if strand == Strands::Reverse {
                exons_.push_front(exon);
            } else {
                exons_.push_back(exon);
            }
        }

        // Step 2. Assign proper exon numbers
        let mut exons: Vec<TranscriptModelExon> = Vec::new();
        let mut exon_number: u16 = 1;
        for exon_ in exons_ {
            let exon: TranscriptModelExon = TranscriptModelExon::new(
                exon_.chromosome_id,
                exon_.start,
                exon_.end,
                exon_number,
                exon_.strand
            );
            exons.push(exon);
            exon_number += 1;
        }
        exons.sort_by_key(|e| e.number);
        exons
    }

    // pub fn identify_overlapping_gene_ids(
    //     exons: &Vec<TranscriptModelExon>,
    //     gene_annotator: &impl GeneAnnotator,
    //     chromosome_names_map: &BiMap<Box<str>,u16>
    // ) -> HashSet<Box<str>> {
    //     let mut gene_ids: HashSet<Box<str>> = HashSet::new();
    //     for exon in exons.iter() {
    //         let chromosome: Box<str> = chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap().to_string().into_boxed_str();
    //         let gene_ids_: Vec<Box<str>> = gene_annotator.get_gene_ids_overlapping_region(&*chromosome, exon.start, exon.end);
    //         for gene_id in gene_ids_ {
    //             let gene: &Gene = gene_annotator.get_gene(&*gene_id).unwrap();
    //             if gene.strand.as_str() == exon.strand.as_str() {
    //                 gene_ids.insert(gene_id);
    //             }
    //         }
    //     }
    //     gene_ids
    // }

    pub fn identify_splice_junctions(
        &self,
        chromosome_names_map: &BiMap<Box<str>,u16>,
        reference_genome_fasta_file: &str,
        min_mapping_quality: usize,
    ) -> Vec<TranscriptModelSpliceJunction> {
        // Step 1. Read the FASTA file
        let mut fasta_reader = Builder::default().build_from_path(reference_genome_fasta_file).unwrap();

        // Step 2. Identify splice junctions
        let mut splice_junctions_: VecDeque<TranscriptModelSpliceJunction> = VecDeque::new();
        for alignment_record in self.alignment_records.iter() {
            // Check if the mapping quality meets the minimum mapping quality
            if min_mapping_quality > alignment_record.record.mapping_quality().unwrap().get() as usize {
                continue;
            }

            // Get the chromosome name
            let chromosome_id: u16 = alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
            let chromosome_name: &str = chromosome_names_map
                .get_by_right(&chromosome_id)
                .unwrap();

            // Get the alignment flag
            let mut strand: Strands = Strands::Forward;
            for flag in alignment_record.record.flags() {
                if flag == Flags::REVERSE_COMPLEMENTED {
                    strand = Strands::Reverse;
                    break;
                }
            }

            // Get the CS tag
            let cs_tag: String;
            if let Some(value) = get_tag_value(&alignment_record.record, "cs") {
                cs_tag = value.to_string();
            } else {
                panic!("Could not find the CS tag.");
            }

            // Get the start reference position
            // curr_reference_pos always points at the end of the last variant
            let mut curr_reference_pos: isize = match alignment_record.record.alignment_start().unwrap() {
                Ok(s) => {
                    s.get() as isize
                },
                Err(e) => {
                    panic!("Could not fetch the start position");
                },
            };
            curr_reference_pos -= 1;

            // Identify splicing in the CS tag
            let re = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
            for value in re.captures_iter(&cs_tag) {
                if value[0].contains(":") {
                    // Increment current reference position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut chars = value[0].chars();
                    chars.next();

                    // Increment current position by the number of matched nucleotides
                    let num_matched_nucleotides: isize = chars.as_str().parse::<isize>().unwrap();
                    curr_reference_pos += num_matched_nucleotides - 1;      // minus 1 is necessary here
                } else if value[0].contains("*") {
                    // Increment current reference position
                    curr_reference_pos += 1;
                } else if value[0].contains("+") {
                    // Do nothing
                } else if value[0].contains("-") {
                    // Increment current position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Reference and alternate alleles
                    let sequence: String = value_chars.as_str().to_string().to_uppercase();
                    let variant_size: usize = sequence.chars().count();

                    // Increment current reference position
                    curr_reference_pos += (variant_size as isize) - 1;
                } else if value[0].contains("~") {
                    // Increment current reference position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Get splicing size
                    let splicing_str: String = value_chars.as_str().to_string();
                    let re_splicing = Regex::new(r"\d+").unwrap(); // r"\d+" matches one or more digits
                    let splicing_size: usize = match re_splicing.find(&splicing_str) {
                        Some(matched) => {
                            let num_str = &splicing_str[matched.start()..matched.end()];
                            let size: usize = match num_str.parse::<usize>() {
                                Ok(num) => num,
                                Err(e) => {
                                    panic!("Failed to convert the splicing size number to usize");
                                },
                            };
                            size
                        },
                        None => {
                            panic!("No numerical value was found in the splicing string: {}", splicing_str);
                        }
                    };

                    // Get the splicing sequence
                    let reference_start: usize =  curr_reference_pos as usize;
                    let reference_end: usize = (curr_reference_pos as usize) + splicing_size - 1;
                    let position_start = Position::try_from(reference_start).unwrap();
                    let position_end = Position::try_from(reference_end).unwrap();
                    let region = Region::new(chromosome_name, position_start..=position_end);
                    let ref_record = fasta_reader.query(&region).unwrap();
                    let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
                    let spliced_sequence_str: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");

                    // Record a splice junction
                    if spliced_sequence_str.len() < 4 {
                        let splice_junction: TranscriptModelSpliceJunction = TranscriptModelSpliceJunction::new(
                            chromosome_id,
                            reference_start as u32,
                            reference_end as u32,
                            0,
                            "".into(),
                            "".into(),
                            strand.clone()
                        );
                        if strand == Strands::Reverse {
                            splice_junctions_.push_front(splice_junction);
                        } else {
                            splice_junctions_.push_back(splice_junction);
                        }
                    } else {
                        let mut splice_signal_start: Box<str> = spliced_sequence_str[..2].to_string().into_boxed_str();
                        let mut splice_signal_end: Box<str> = spliced_sequence_str[spliced_sequence_str.len()-2..].to_string().into_boxed_str();
                        if strand == Strands::Reverse {
                            splice_signal_start = reverse_complement(&*splice_signal_start);
                            splice_signal_end = reverse_complement(&*splice_signal_end);
                        }
                        let splice_junction: TranscriptModelSpliceJunction = TranscriptModelSpliceJunction::new(
                            chromosome_id,
                            reference_start as u32,
                            reference_end as u32,
                            0,
                            &*splice_signal_start,
                            &*splice_signal_end,
                            strand.clone()
                        );
                        if strand == Strands::Reverse {
                            splice_junctions_.push_front(splice_junction);
                        } else {
                            splice_junctions_.push_back(splice_junction);
                        }
                    }

                    // Increment current reference position
                    curr_reference_pos += (splicing_size as isize) - 1;
                } else {
                    panic!("Unknown delimiter for CS tag: {}", value[0].to_string());
                }
            }
        }

        // Step 2. Assign proper splice junction numbers
        let mut splice_junctions: Vec<TranscriptModelSpliceJunction> = Vec::new();
        let mut splice_junction_number: u16 = 1;
        for splice_junction_ in splice_junctions_ {
            let splice_junction: TranscriptModelSpliceJunction = TranscriptModelSpliceJunction::new(
                splice_junction_.chromosome_id,
                splice_junction_.start,
                splice_junction_.end,
                splice_junction_number,
                &*splice_junction_.splice_signal_start,
                &*splice_junction_.splice_signal_end,
                splice_junction_.strand
            );
            splice_junctions.push(splice_junction);
            splice_junction_number += 1;
        }

        splice_junctions
    }

    pub fn identify_splice_variant_records(
        &self, 
        matched_reference_transcripts: &Vec<ReferenceTranscriptMatch>,
        chromosome_names_map: &BiMap<Box<str>,u16>,
        min_mapping_quality: usize
    ) -> Vec<VariantRecord> {
        // Step 1. Identify exons
        let exons: Vec<TranscriptModelExon> = self.identify_exons(min_mapping_quality);
        
        // Step 2. Identify splice variant records
        let mut splice_variant_records: Vec<VariantRecord> = Vec::new();
        if matched_reference_transcripts.is_empty() {
            // All the exons will be identified as cryptic exons
            for exon in exons.iter() {
                let sequence_operation: SequenceOperation = SequenceOperation::new(
                    exon.chromosome_id,
                    exon.start,
                    exon.strand.clone(),
                    SequenceOperationTypes::Read.clone(),
                    exon.chromosome_id,
                    exon.end,
                    exon.strand.clone(),
                    SequenceOperationTypes::Read.clone(),
                    "".into(),
                    SequenceOperationVariantTypes::CrypticExon.clone()
                );
                let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                splice_variant_records.push(variant_record);
            }
        } else {
            let mut reference_transcripts: Vec<Transcript> = Vec::new();
            for matched_reference_transcript in matched_reference_transcripts.iter() {
                reference_transcripts.push(matched_reference_transcript.reference_transcript.clone());
            }

            // Identify fusion
            let mut curr_gene_id: Box<str> = "".into();
            let mut prev_exon: &TranscriptModelExon = exons.first().unwrap();
            for (index,exon) in exons.iter().enumerate() {
                // Found the overlapping gene ID
                let exon_start: isize = exon.start as isize;
                let exon_end: isize = exon.end as isize;
                let mut overlapping_gene_id: Box<str> = "".into();
                for reference_transcript in reference_transcripts.iter() {
                    let reference_transcript_start: isize = reference_transcript.start as isize;
                    let reference_transcript_end: isize = reference_transcript.end as isize;
                    let reference_chromosome_id: u16 = *chromosome_names_map.get_by_left(&*reference_transcript.chromosome).unwrap();
                    if reference_chromosome_id == exon.chromosome_id &&
                        exon.strand.as_str() == reference_transcript.strand.as_str() {
                        if overlaps(exon_start, exon_end, reference_transcript_start, reference_transcript_end) {
                            overlapping_gene_id = reference_transcript.gene_id.clone();
                            break;
                        }
                    }
                }
                if index == 0 {
                    if overlapping_gene_id != "".into() {
                        curr_gene_id = overlapping_gene_id;
                    }
                } else {
                    if overlapping_gene_id != curr_gene_id && overlapping_gene_id != "".into() {
                        if prev_exon.strand.as_str() == Strands::Forward.as_str() {
                            let sequence_operation: SequenceOperation= SequenceOperation::new(
                                prev_exon.chromosome_id,
                                prev_exon.end,
                                prev_exon.strand.clone(),
                                SequenceOperationTypes::Downstream.clone(),
                                exon.chromosome_id,
                                exon.start,
                                exon.strand.clone(),
                                SequenceOperationTypes::Upstream.clone(),
                                "".into(),
                                SequenceOperationVariantTypes::FusionGene.clone()
                            );
                            let variant_record: VariantRecord = VariantRecord::new(
                                self.read_id, sequence_operation
                            );
                            splice_variant_records.push(variant_record);
                        } else {
                            let sequence_operation: SequenceOperation = SequenceOperation::new(
                                prev_exon.chromosome_id,
                                prev_exon.start,
                                prev_exon.strand.clone(),
                                SequenceOperationTypes::Upstream.clone(),
                                exon.chromosome_id,
                                exon.end,
                                exon.strand.clone(),
                                SequenceOperationTypes::Downstream.clone(),
                                "".into(),
                                SequenceOperationVariantTypes::FusionGene.clone()
                            );
                            let variant_record: VariantRecord = VariantRecord::new(
                                self.read_id, sequence_operation
                            );
                            splice_variant_records.push(variant_record);
                        }
                        curr_gene_id = overlapping_gene_id;
                    }
                }
                prev_exon = exon;
            }

            // Identify reference exon skipping
            for reference_transcript in reference_transcripts.iter() {
                for reference_exon in reference_transcript.exons.values() {
                    let mut reference_exon_overlaps: bool = false;
                    let reference_exon_start: isize = reference_exon.start as isize;
                    let reference_exon_end: isize = reference_exon.end as isize;
                    let reference_chromosome_id: u16 = *chromosome_names_map.get_by_left(&*reference_exon.chromosome).unwrap();
                    let reference_strand: Strands = if reference_exon.strand.as_str() == Strands::Forward.as_str() {
                        Strands::Forward
                    } else {
                        Strands::Reverse
                    };
                    for exon in exons.iter() {
                        if reference_chromosome_id == exon.chromosome_id {
                            let exon_start: isize = exon.start as isize;
                            let exon_end: isize = exon.end as isize;
                            if overlaps(exon_start, exon_end, reference_exon_start, reference_exon_end) {
                                reference_exon_overlaps = true;
                                break;
                            }
                        }
                    }

                    // Check if the reference exon is skipped because of a fusion gene
                    for variant_record in splice_variant_records.iter() {
                        if variant_record.get_variant_type() == SequenceOperationVariantTypes::FusionGene {
                            let position_1: isize = variant_record.get_position_1() as isize;
                            let position_2: isize = variant_record.get_position_2() as isize;
                            if reference_chromosome_id == variant_record.get_chromosome_1() && reference_chromosome_id == variant_record.get_chromosome_2() {
                                if overlaps(position_1, position_2, reference_exon_start, reference_exon_end) {
                                    reference_exon_overlaps = true;
                                    break;
                                }
                            }
                        }
                    }

                    if reference_exon_overlaps == false {
                        let strand = if reference_exon.strand.as_str() == Strands::Forward.as_str() {
                            Strands::Forward
                        } else {
                            Strands::Reverse
                        };
                        let sequence_operation: SequenceOperation = SequenceOperation::new(
                            reference_chromosome_id,
                            reference_exon.start,
                            strand.clone(),
                            SequenceOperationTypes::Skip.clone(),
                            reference_chromosome_id,
                            reference_exon.end,
                            strand.clone(),
                            SequenceOperationTypes::Skip.clone(),
                            "".into(),
                            SequenceOperationVariantTypes::ExonSkipping.clone()
                        );
                        let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                        splice_variant_records.push(variant_record);
                    }
                }
            }

            // Identify cryptic exons
            for exon in exons.iter() {
                let mut is_cryptic: bool = true;
                let exon_start: isize = exon.start as isize;
                let exon_end: isize = exon.end as isize;
                let chromosome_name: Box<str> = chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap().to_string().into_boxed_str();
                for reference_transcript in reference_transcripts.iter() {
                    for reference_exon in reference_transcript.exons.values() {
                        if chromosome_name == reference_exon.chromosome {
                            let reference_exon_start: isize = reference_exon.start as isize;
                            let reference_exon_end: isize = reference_exon.end as isize;
                            if overlaps(exon_start, exon_end, reference_exon_start, reference_exon_end) {
                                is_cryptic = false;
                                break;
                            }
                        }
                    }
                    if is_cryptic == false {
                        break;
                    }
                }
                if is_cryptic {
                    let sequence_operation: SequenceOperation = SequenceOperation::new(
                        exon.chromosome_id,
                        exon.start,
                        exon.strand.clone(),
                        SequenceOperationTypes::Read.clone(),
                        exon.chromosome_id,
                        exon.end,
                        exon.strand.clone(),
                        SequenceOperationTypes::Read.clone(),
                        "".into(),
                        SequenceOperationVariantTypes::CrypticExon.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                    splice_variant_records.push(variant_record);
                }
            }

            // Identify intron retention
            for exon in exons.iter() {
                let mut overlaps_exon: bool = false;
                let mut overlaps_intron: bool = false;
                let mut intron_retention_start: u32 = 0;
                let mut intron_retention_end: u32 = 0;
                let exon_start: isize = exon.start as isize;
                let exon_end: isize = exon.end as isize;
                for reference_transcript in reference_transcripts.iter() {
                    for intron in reference_transcript.get_introns().iter() {
                        let intron_start: isize = intron.start as isize;
                        let intron_end: isize = intron.end as isize;
                        let intron_chromosome_id: u16 = *chromosome_names_map.get_by_left(&intron.chromosome).unwrap();
                        if intron_chromosome_id == exon.chromosome_id &&
                            overlaps(exon_start, exon_end, intron_start, intron_end) {
                            let (start, end) = find_overlap((exon_start, exon_end), (intron_start, intron_end)).unwrap();
                            overlaps_intron = true;
                            intron_retention_start = start as u32;
                            intron_retention_end = end as u32;
                            break;
                        }
                    }
                    for reference_exon in reference_transcript.exons.values() {
                        let reference_exon_start: isize = reference_exon.start as isize;
                        let reference_exon_end: isize = reference_exon.end as isize;
                        let exon_chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_exon.chromosome).unwrap();
                        if exon_chromosome_id == exon.chromosome_id &&
                            overlaps(exon_start, exon_end, reference_exon_start, reference_exon_end) {
                            overlaps_exon = true;
                            break;
                        }
                    }
                }
                if overlaps_intron && overlaps_exon {
                    let sequence_operation: SequenceOperation = SequenceOperation::new(
                        exon.chromosome_id,
                        intron_retention_start,
                        exon.strand.clone(),
                        SequenceOperationTypes::Read.clone(),
                        exon.chromosome_id,
                        intron_retention_end,
                        exon.strand.clone(),
                        SequenceOperationTypes::Read.clone(),
                        "".into(),
                        SequenceOperationVariantTypes::IntronRetention.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                    splice_variant_records.push(variant_record);
                }
            }

            // Identify alternative splice sites
            for reference_transcript in reference_transcripts.iter() {
                for reference_exon in reference_transcript.get_sorted_exons() {
                    let reference_exon_start: isize = reference_exon.start as isize;
                    let reference_exon_end: isize = reference_exon.end as isize;
                    let reference_chromosome_id: u16 = *chromosome_names_map.get_by_left(&reference_exon.chromosome).unwrap();
                    for exon in exons.iter() {
                        let exon_start: isize = exon.start as isize;
                        let exon_end: isize = exon.end as isize;
                        if reference_chromosome_id != exon.chromosome_id {
                            continue;
                        }
                        if overlaps(exon_start, exon_end, reference_exon_start, reference_exon_end) == false {
                            continue;
                        }
                        let (start,end) = find_overlap((exon_start, exon_end), (reference_exon_start, reference_exon_end)).unwrap();

                        // Check if the alternative splice site is a fusion gene breakpoint
                        let mut is_fusion_breakpoint: bool = false;
                        for variant_record in splice_variant_records.iter() {
                            if variant_record.get_variant_type() == SequenceOperationVariantTypes::FusionGene {
                                if exon.chromosome_id == variant_record.get_chromosome_1() && exon.chromosome_id == variant_record.get_chromosome_2() {
                                    let position_1: isize = variant_record.get_position_1() as isize;
                                    let position_2: isize = variant_record.get_position_2() as isize;
                                    if start >= position_1 - 1 && start <= position_1 + 1 {
                                        is_fusion_breakpoint = true;
                                    }
                                    if start >= position_2 - 1 && start <= position_2 + 1 {
                                        is_fusion_breakpoint = true;
                                    }
                                    if end >= position_1 - 1 && end <= position_1 + 1 {
                                        is_fusion_breakpoint = true;
                                    }
                                    if end >= position_2 - 1 && end <= position_2 + 1 {
                                        is_fusion_breakpoint = true;
                                    }
                                }
                            }
                        }
                        if start != reference_exon_start && is_fusion_breakpoint == false {
                            // Check if 5' or 3' splice site
                            if reference_exon.strand.as_str() == Strands::Forward.as_str() {
                                let sequence_operation: SequenceOperation = SequenceOperation::new(
                                    exon.chromosome_id,
                                    exon.start - 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    exon.chromosome_id,
                                    exon.start - 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    "".into(),
                                    SequenceOperationVariantTypes::Alternative3PrimeSpliceSite.clone()
                                );
                                let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                                splice_variant_records.push(variant_record);
                            } else {
                                let sequence_operation: SequenceOperation = SequenceOperation::new(
                                    exon.chromosome_id,
                                    exon.start - 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    exon.chromosome_id,
                                    exon.start - 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    "".into(),
                                    SequenceOperationVariantTypes::Alternative5PrimeSpliceSite.clone()
                                );
                                let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                                splice_variant_records.push(variant_record);
                            }
                        }
                        if end != reference_exon_end && is_fusion_breakpoint == false {
                            // Check if 5' or 3' splice site
                            if reference_exon.strand.as_str() == Strands::Forward.as_str() {
                                let sequence_operation: SequenceOperation = SequenceOperation::new(
                                    exon.chromosome_id,
                                    exon.end + 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    exon.chromosome_id,
                                    exon.end + 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    "".into(),
                                    SequenceOperationVariantTypes::Alternative5PrimeSpliceSite.clone()
                                );
                                let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                                splice_variant_records.push(variant_record);
                            } else {
                                let sequence_operation: SequenceOperation = SequenceOperation::new(
                                    exon.chromosome_id,
                                    exon.end + 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    exon.chromosome_id,
                                    exon.end + 1,
                                    exon.strand.clone(),
                                    SequenceOperationTypes::Mark.clone(),
                                    "".into(),
                                    SequenceOperationVariantTypes::Alternative3PrimeSpliceSite.clone()
                                );
                                let variant_record: VariantRecord = VariantRecord::new(self.read_id, sequence_operation);
                                splice_variant_records.push(variant_record);
                            }
                        }
                    }
                }
            }
        }

        splice_variant_records
    }

    pub fn identify_sequence_variant_records(
        &self,
        min_mapping_quality: usize,
        min_average_base_quality: f32
    ) -> Vec<VariantRecord> {
        let mut variant_records = self.identify_sequence_variant_records_in_cs_tag(
            min_mapping_quality,
            min_average_base_quality
        );
        variant_records.extend(
            self.identify_sequence_variant_records_in_softclipping(
                min_mapping_quality,
                min_average_base_quality
            )
        );
        variant_records
    }

    /// Identify intra-chromosomal SNVs, insertions, splicing, and deletions using CS tag.
    pub fn identify_sequence_variant_records_in_cs_tag(
        &self,
        min_mapping_quality: usize,
        min_average_base_quality: f32
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        for alignment_record in self.alignment_records.iter() {
            // Check if the mapping quality meets the minimum mapping quality
            if min_mapping_quality > alignment_record.record.mapping_quality().unwrap().get() as usize {
                continue;
            }

            // Get the alignment base quality scores
            let base_quality_scores: Vec<u8> = get_alignment_base_quality_scores(&alignment_record.record);
            let mut curr_base_quality_pos: isize = -1;

            // Get the alignment flag
            let mut strand: Strands = Strands::Forward;
            for flag in alignment_record.record.flags() {
                if flag == Flags::REVERSE_COMPLEMENTED {
                    strand = Strands::Reverse;
                    break;
                }
            }

            // Get the CS tag
            let cs_tag: String;
            if let Some(value) = get_tag_value(&alignment_record.record, "cs") {
                cs_tag = value.to_string();
            } else {
                panic!("Could not find the CS tag.");
            }

            // Get the start reference position
            // curr_reference_pos always points at the end of the last variant
            let mut curr_reference_pos: isize = match alignment_record.record.alignment_start().unwrap() {
                Ok(s) => {
                    s.get() as isize
                },
                Err(e) => {
                    panic!("Could not fetch the start position");
                },
            };
            curr_reference_pos -= 1;

            // Identify SNVs, insertions, deletions, and splicing in the CS tag
            // let chromosome_name: &str = chromosome_names_map
            //     .get_by_right(&alignment_record.record.reference_sequence_id().unwrap().unwrap())
            //     .unwrap();
            let chromosome_id: u16 = alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
            let re = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
            for value in re.captures_iter(&cs_tag) {
                if value[0].contains(":") {
                    // Increment current reference position
                    curr_reference_pos += 1;
                    curr_base_quality_pos += 1;

                    // Remove the first character
                    let mut chars = value[0].chars();
                    chars.next();

                    // Increment current position by the number of matched nucleotides
                    let num_matched_nucleotides: isize = chars.as_str().parse::<isize>().unwrap();
                    curr_reference_pos += num_matched_nucleotides - 1;      // minus 1 is necessary here
                    curr_base_quality_pos += num_matched_nucleotides - 1;   // minus 1 is necessary here
                } else if value[0].contains("*") {
                    // Increment current reference position
                    curr_reference_pos += 1;
                    curr_base_quality_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Reference and alternate alleles
                    let alleles: String = value_chars.as_str().to_string();
                    let reference_allele: String = alleles.chars().nth(0).unwrap().to_string().to_uppercase();
                    let alternate_allele: String = alleles.chars().nth(1).unwrap().to_string().to_uppercase();

                    // Check base quality score
                    let sequence_quality_scores: Vec<u8> = vec![base_quality_scores[curr_base_quality_pos as usize]];
                    if min_average_base_quality > sequence_quality_scores[0] as f32 {
                        continue;
                    }

                    // Exclude SNVs at the first or last base of the read
                    if curr_base_quality_pos == 0 || curr_base_quality_pos == self.get_read_length() as isize - 1 {
                        continue;
                    }

                    // Record a single-nucleotide variant
                    let graph_operation: SequenceOperation = SequenceOperation::new(
                        chromosome_id,
                        (curr_reference_pos - 1) as u32,
                        strand.clone(),
                        SequenceOperationTypes::Downstream,
                        chromosome_id,
                        (curr_reference_pos + 1) as u32,
                        strand.clone(),
                        SequenceOperationTypes::Upstream,
                        alternate_allele.into_boxed_str(),
                        SequenceOperationVariantTypes::SingleNucleotideVariant.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(
                        self.read_id,
                        graph_operation
                    );
                    variant_records.push(variant_record);
                } else if value[0].contains("+") {
                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Reference and alternate alleles
                    let alternate_allele: String = value_chars.as_str().to_string().to_uppercase();
                    let variant_size: usize = alternate_allele.chars().count();

                    // Base quality score
                    curr_base_quality_pos += 1;
                    let curr_base_quality_start: usize = curr_base_quality_pos as usize;
                    let curr_base_quality_end: usize = curr_base_quality_start + variant_size - 1;
                    let sequence_quality_scores: Vec<u8> = base_quality_scores[curr_base_quality_start..=curr_base_quality_end].to_vec();
                    curr_base_quality_pos += (variant_size as isize) - 1;
                    if min_average_base_quality > calculate_average_base_quality_score(&sequence_quality_scores) {
                        continue;
                    }

                    // Record an insertion
                    let graph_operation: SequenceOperation = SequenceOperation::new(
                        chromosome_id,
                        curr_reference_pos as u32,
                        strand.clone(),
                        SequenceOperationTypes::Downstream,
                        chromosome_id,
                        (curr_reference_pos + 1) as u32,
                        strand.clone(),
                        SequenceOperationTypes::Upstream,
                        alternate_allele.into_boxed_str(),
                        SequenceOperationVariantTypes::Insertion.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(
                        self.read_id,
                        graph_operation
                    );
                    variant_records.push(variant_record);
                } else if value[0].contains("-") {
                    // Increment current position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Reference and alternate alleles
                    let reference_allele: String = value_chars.as_str().to_string().to_uppercase();
                    let alternate_allele: String = String::from("");
                    let sequence: String = value_chars.as_str().to_string().to_uppercase();
                    let variant_size: usize = sequence.chars().count() as usize;

                    // Get the deletion sequence
                    let reference_start: usize = curr_reference_pos as usize;
                    let reference_end: usize = (curr_reference_pos as usize) + variant_size - 1;
                    // let position_start: Position = Position::try_from(reference_start).unwrap();
                    // let position_end: Position = Position::try_from(reference_end).unwrap();
                    // let region: Region = Region::new(chromosome_name, position_start..=position_end);
                    // let ref_record: fasta::Record = reader.query(&region).unwrap();
                    // let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
                    // let deletion_sequence_str: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");

                    // Record a deletion
                    let graph_operation: SequenceOperation = SequenceOperation::new(
                        chromosome_id,
                        (reference_start - 1) as u32,
                        strand.clone(),
                        SequenceOperationTypes::Downstream,
                        chromosome_id,
                        (reference_end + 1) as u32,
                        strand.clone(),
                        SequenceOperationTypes::Upstream,
                        "".to_string().into_boxed_str(),
                        SequenceOperationVariantTypes::Deletion.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(
                        self.read_id,
                        graph_operation
                    );
                    variant_records.push(variant_record);

                    // Increment current reference position
                    curr_reference_pos += (variant_size as isize) - 1;
                } else if value[0].contains("~") {
                    // Increment current reference position
                    curr_reference_pos += 1;

                    // Remove the first character
                    let mut value_chars = value[0].chars();
                    value_chars.next();

                    // Get splicing size
                    let splicing_str: String = value_chars.as_str().to_string();
                    let re_splicing = Regex::new(r"\d+").unwrap(); // r"\d+" matches one or more digits
                    let splicing_size: usize = match re_splicing.find(&splicing_str) {
                        Some(matched) => {
                            let num_str = &splicing_str[matched.start()..matched.end()];
                            let size: usize = match num_str.parse::<usize>() {
                                Ok(num) => num,
                                Err(e) => {
                                    panic!("Failed to convert the splicing size number to usize");
                                },
                            };
                            size
                        },
                        None => {
                            panic!("No numerical value was found in the splicing string: {}", splicing_str);
                        }
                    };

                    // Get the splicing sequence
                    let reference_start: usize =  curr_reference_pos as usize;
                    let reference_end: usize = (curr_reference_pos as usize) + splicing_size - 1;
                    // let position_start = Position::try_from(start).unwrap();
                    // let position_end = Position::try_from(end).unwrap();
                    // let region = Region::new(chromosome, position_start..=position_end);
                    // let ref_record = reader.query(&region).unwrap();
                    // let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
                    // let deletion_sequence_str: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");

                    // Record a splicing
                    // let graph_operation: GraphOperation = GraphOperation::new(
                    //     chromosome_id,
                    //     (reference_start - 1) as u32,
                    //     strand.clone(),
                    //     GraphOperationOrientations::DOWNSTREAM,
                    //     chromosome_id,
                    //     (reference_end + 1) as u32,
                    //     strand.clone(),
                    //     GraphOperationOrientations::UPSTREAM,
                    //     "".to_string().into_boxed_str()
                    // );
                    // let variant_record: VariantRecord = VariantRecord::new(
                    //     self.read_id,
                    //     graph_operation
                    // );
                    // variant_records.push(variant_record);

                    // Increment current reference position
                    curr_reference_pos += (splicing_size as isize) - 1;
                } else {
                    panic!("Unknown delimiter for CS tag: {}", value[0].to_string());
                }
            }
        }

        variant_records
    }

    /// Identify variant records using softclipped bases.
    ///
    /// # Parameters
    ///
    /// * `records` is a vector of all BAM records of the same read ID.
    /// * `chromosomes` is a HashMap where `key` is chromosome ID and `value` is chromosome name.
    pub fn identify_sequence_variant_records_in_softclipping(
        &self,
        min_mapping_quality: usize,
        min_average_base_quality: f32
    ) -> Vec<VariantRecord> {
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        let mut prev_alignment_record: &AlignmentRecord = &self.alignment_records[0];
        for (i,curr_alignment_record) in self.alignment_records.iter().enumerate() {
            // Identify soft-clipped insertion in the first alignment
            if (i == 0) && (curr_alignment_record.read_start != 0) {
                let chromosome_id: u16 = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                let read_sequence: Box<str> = get_read_sequence(&curr_alignment_record.record);
                let base_quality_scores: Vec<u8> = get_base_quality_scores(&curr_alignment_record.record);
                let position_1: usize;
                let position_2: usize;
                let operation_1: SequenceOperationTypes;
                let operation_2: SequenceOperationTypes;
                let insertion: String;
                let sequence_quality_scores: Vec<u8>;

                if curr_alignment_record.record.flags().is_reverse_complemented() {
                    let right_softclipping: (bool,usize) = get_right_softclipping(&curr_alignment_record.record);
                    assert!(
                        right_softclipping.0 && right_softclipping.1 == curr_alignment_record.read_start,
                        "The 3' end of the first alignment (read ID: {}) is softclipped so the alignment's read start position is expected to be the same as the number of softclipped bases.", self.read_id
                    );
                    position_1 = get_alignment_end_position(&curr_alignment_record.record);
                    position_2 = get_alignment_end_position(&curr_alignment_record.record) + 1;
                    operation_1 = SequenceOperationTypes::Downstream;
                    operation_2 = SequenceOperationTypes::Upstream;
                    insertion = read_sequence[self.get_read_length() - curr_alignment_record.read_start..self.get_read_length()].to_string();
                    sequence_quality_scores = base_quality_scores[self.get_read_length() - curr_alignment_record.read_start..self.get_read_length()].to_vec();
                } else {
                    let left_softclipping: (bool,usize) = get_left_softclipping(&curr_alignment_record.record);
                    assert!(
                        left_softclipping.0 && left_softclipping.1 == curr_alignment_record.read_start,
                        "The 5' end of the first alignment (read ID: {}) is softclipped so the alignment's read start position is expected to be the same as the number of softclipped bases.", self.read_id
                    );
                    position_1 = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    position_2 = get_alignment_start_position(&curr_alignment_record.record);
                    operation_1 = SequenceOperationTypes::Downstream;
                    operation_2 = SequenceOperationTypes::Upstream;
                    insertion = read_sequence[0..curr_alignment_record.read_start].to_string();
                    sequence_quality_scores = base_quality_scores[0..curr_alignment_record.read_start].to_vec();
                }

                let strand: Strands = if is_aligned_to_reverse_strand(&curr_alignment_record.record) {
                    Strands::Reverse
                } else {
                    Strands::Forward
                };

                // Exclude the following soft-clipping cases:
                // - Does not meet the minimum mapping quality
                // - Does not meet the minimum average base quality
                // - Single base soft-clipping
                // - 2 or 3-base soft-clipping with AA, CC, GG, or TT
                let mut include: bool = true;
                if min_mapping_quality > curr_alignment_record.record.mapping_quality().unwrap().get() as usize {
                    include = false;
                }
                if min_average_base_quality > calculate_average_base_quality_score(&sequence_quality_scores) {
                    include = false;
                }
                if insertion.len() == 1 {
                    include = false;
                }
                if insertion.len() <= 3 {
                    if insertion.to_uppercase().contains("AA") ||
                        insertion.to_uppercase().contains("CC") ||
                        insertion.to_uppercase().contains("GG") ||
                        insertion.to_uppercase().contains("TT") {
                        include = false;
                    }
                }
                if include {
                    let graph_operation: SequenceOperation = SequenceOperation::new(
                        chromosome_id,
                        position_1 as u32,
                        strand.clone(),
                        operation_1.clone(),
                        chromosome_id,
                        position_2 as u32,
                        strand.clone(),
                        operation_1.clone(),
                        insertion.into_boxed_str(),
                        SequenceOperationVariantTypes::Insertion.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(
                        self.read_id,
                        graph_operation
                    );
                    variant_records.push(variant_record);
                }
            }

            // Identify soft-clipped insertion in the last alignment
            if (i == self.alignment_records.len() - 1) && ((curr_alignment_record.read_end) != self.get_read_length() - 1) {
                let chromosome_id: u16 = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                let read_sequence: Box<str> = get_read_sequence(&curr_alignment_record.record);
                let base_quality_scores: Vec<u8> = get_base_quality_scores(&curr_alignment_record.record);
                let position_1: usize;
                let position_2: usize;
                let operation_1: SequenceOperationTypes;
                let operation_2: SequenceOperationTypes;
                let insertion: String;
                let sequence_quality_scores: Vec<u8>;
                if curr_alignment_record.record.flags().is_reverse_complemented() {
                    let left_softclipping: (bool,usize) = get_left_softclipping(&curr_alignment_record.record);
                    assert!(
                        left_softclipping.0 && left_softclipping.1 == (self.get_read_length() - curr_alignment_record.read_end - 1),
                        "The 5' end of the last alignment (read ID: {}) is softclipped so (read length - alignment's last read position - 1) is expected to be the same as the number of softclipped bases.", self.read_id
                    );
                    position_1 = get_alignment_start_position(&curr_alignment_record.record) - 1;
                    position_2 = get_alignment_start_position(&curr_alignment_record.record);
                    operation_1 = SequenceOperationTypes::Downstream;
                    operation_2 = SequenceOperationTypes::Upstream;
                    insertion = read_sequence[0..(self.get_read_length() - curr_alignment_record.read_end - 1)].to_string();
                    sequence_quality_scores = base_quality_scores[0..(self.get_read_length() - curr_alignment_record.read_end - 1)].to_vec();
                } else {
                    let right_softclipping: (bool,usize) = get_right_softclipping(&curr_alignment_record.record);
                    assert!(
                        right_softclipping.0 && right_softclipping.1 == (self.get_read_length() - curr_alignment_record.read_end - 1),
                        "The 3' end of the last alignment (read ID: {}) is softclipped so (read length - alignment's last read position - 1) is expected to be the same as the number of softclipped bases.", self.read_id
                    );
                    position_1 = get_alignment_end_position(&curr_alignment_record.record);
                    position_2 = get_alignment_end_position(&curr_alignment_record.record) + 1;
                    operation_1 = SequenceOperationTypes::Downstream;
                    operation_2 = SequenceOperationTypes::Upstream;
                    insertion = read_sequence[curr_alignment_record.read_end+1..self.get_read_length()].to_string();
                    sequence_quality_scores = base_quality_scores[(curr_alignment_record.read_end + 1)..self.get_read_length()].to_vec();
                }
                let strand: Strands = if is_aligned_to_reverse_strand(&curr_alignment_record.record) {
                    Strands::Reverse
                } else {
                    Strands::Forward
                };

                // Exclude the following soft-clipping cases:
                // - Does not meet the minimum mapping quality
                // - Does not meet the minimum average base quality
                // - Single base soft-clipping
                // - 2 or 3-base soft-clipping with AA, CC, GG, or TT
                let mut include: bool = true;
                if min_mapping_quality > curr_alignment_record.record.mapping_quality().unwrap().get() as usize {
                    include = false;
                }
                if min_average_base_quality > calculate_average_base_quality_score(&sequence_quality_scores) {
                    include = false;
                }
                if insertion.len() == 1 {
                    include = false;
                }
                if insertion.len() <= 3 {
                    if insertion.to_uppercase().contains("AA") ||
                        insertion.to_uppercase().contains("CC") ||
                        insertion.to_uppercase().contains("GG") ||
                        insertion.to_uppercase().contains("TT") {
                        include = false;
                    }
                }
                if include {
                    let graph_operation: SequenceOperation = SequenceOperation::new(
                        chromosome_id,
                        position_1 as u32,
                        strand.clone(),
                        operation_1.clone(),
                        chromosome_id,
                        position_2 as u32,
                        strand.clone(),
                        operation_2.clone(),
                        insertion.into_boxed_str(),
                        SequenceOperationVariantTypes::Insertion.clone()
                    );
                    let variant_record: VariantRecord = VariantRecord::new(
                        self.read_id,
                        graph_operation
                    );
                    variant_records.push(variant_record);
                }
            }

            // Identify breakpoints (softclipping) between alignments
            let prev_alignment_mapping_quality: usize = prev_alignment_record.record.mapping_quality().unwrap().get() as usize;
            let curr_alignment_mapping_quality: usize = curr_alignment_record.record.mapping_quality().unwrap().get() as usize;
            if i > 0 &&
                prev_alignment_mapping_quality >= min_mapping_quality &&
                curr_alignment_mapping_quality >= min_mapping_quality {
                // Determine the following
                let bnd_chromosome_1_id: u16;
                let bnd_chromosome_2_id: u16;
                let mut bnd_position_1: usize;
                let mut bnd_position_2: usize;
                let bnd_operation_1: SequenceOperationTypes;
                let bnd_operation_2: SequenceOperationTypes;
                let mut insertion: Box<str> = "".to_string().into_boxed_str();
                let mut insertion_sequence_quality_scores: Vec<u8> = vec![];
                if !prev_alignment_record.record.flags().is_reverse_complemented() && curr_alignment_record.record.flags().is_reverse_complemented() {
                    bnd_chromosome_1_id = prev_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_chromosome_2_id = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_position_1 = prev_alignment_record.record.alignment_end().unwrap().unwrap().get();
                    bnd_position_2 = curr_alignment_record.record.alignment_end().unwrap().unwrap().get();
                    bnd_operation_1 = SequenceOperationTypes::Downstream;
                    bnd_operation_2 = SequenceOperationTypes::Downstream;
                } else if prev_alignment_record.record.flags().is_reverse_complemented() && !curr_alignment_record.record.flags().is_reverse_complemented() {
                    bnd_chromosome_1_id = prev_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_chromosome_2_id = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_position_1 = prev_alignment_record.record.alignment_start().unwrap().unwrap().get();
                    bnd_position_2 = curr_alignment_record.record.alignment_start().unwrap().unwrap().get();
                    bnd_operation_1 = SequenceOperationTypes::Upstream;
                    bnd_operation_2 = SequenceOperationTypes::Upstream;
                } else if !prev_alignment_record.record.flags().is_reverse_complemented() && !curr_alignment_record.record.flags().is_reverse_complemented() {
                    bnd_chromosome_1_id = prev_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_chromosome_2_id = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_position_1 = prev_alignment_record.record.alignment_end().unwrap().unwrap().get();
                    bnd_position_2 = curr_alignment_record.record.alignment_start().unwrap().unwrap().get();
                    bnd_operation_1 = SequenceOperationTypes::Downstream;
                    bnd_operation_2 = SequenceOperationTypes::Upstream;
                } else {
                    bnd_chromosome_1_id = prev_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_chromosome_2_id = curr_alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
                    bnd_position_1 = prev_alignment_record.record.alignment_start().unwrap().unwrap().get();
                    bnd_position_2 = curr_alignment_record.record.alignment_end().unwrap().unwrap().get();
                    bnd_operation_1 = SequenceOperationTypes::Upstream;
                    bnd_operation_2 = SequenceOperationTypes::Downstream;
                }

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
                    insertion = self.original_read_sequence[(overlap_start as usize)..=(overlap_end as usize)].to_string().into_boxed_str();
                    insertion_sequence_quality_scores = self.quality_scores[(overlap_start as usize)..=(overlap_end as usize)].to_vec();

                    if prev_alignment_record.reverse_complemented && curr_alignment_record.reverse_complemented {
                        insertion = reverse_complement(&*insertion);
                        insertion_sequence_quality_scores.reverse();
                    }

                    // Update position 1 by the number of overlapping bases
                    if bnd_operation_1 == SequenceOperationTypes::Upstream {
                        bnd_position_1 = bnd_position_1 + insertion.len();
                    } else {
                       bnd_position_1 = bnd_position_1 - insertion.len();
                    }

                    // Update position 2 by the number of overlapping bases
                    if bnd_operation_2 == SequenceOperationTypes::Upstream {
                        bnd_position_2 = bnd_position_2 + insertion.len();
                    } else {
                        bnd_position_2 = bnd_position_2 - insertion.len();
                    }
                } else {
                    // Check if an insertion exists between the breakpoints
                    if prev_alignment_record.read_end+1 != curr_alignment_record.read_start &&
                        prev_alignment_record.read_end < curr_alignment_record.read_start {
                        insertion = self.original_read_sequence[(prev_alignment_record.read_end+1)..curr_alignment_record.read_start].to_string().into_boxed_str();
                        insertion_sequence_quality_scores = self.quality_scores[(prev_alignment_record.read_end+1)..curr_alignment_record.read_start].to_vec();
                        if prev_alignment_record.reverse_complemented && curr_alignment_record.reverse_complemented {
                            insertion = reverse_complement(&*insertion);
                            insertion_sequence_quality_scores.reverse();
                        }
                    }
                }

                // Exclude the insertion if it does not meet the minimum average base quality
                if insertion != "".into() {
                    if min_average_base_quality > calculate_average_base_quality_score(&insertion_sequence_quality_scores) {
                        insertion = "".to_string().into_boxed_str();
                    }
                }

                // Determine the strand
                let bnd_strand_1: Strands = if is_aligned_to_reverse_strand(&prev_alignment_record.record) {
                    Strands::Reverse
                } else {
                    Strands::Forward
                };
                let bnd_strand_2: Strands = if is_aligned_to_reverse_strand(&curr_alignment_record.record) {
                    Strands::Reverse
                } else {
                    Strands::Forward
                };

                // Enforce an order for chromosomes 1 and 2 as well as positions 1 and 2
                let (chromosome_1_id,position_1,strand_1,orientation_1,chromosome_2_id,position_2,strand_2,orientation_2) =
                    if bnd_chromosome_1_id < bnd_chromosome_2_id ||
                        (bnd_chromosome_1_id == bnd_chromosome_2_id && bnd_position_1 < bnd_position_2) {
                        (
                            bnd_chromosome_1_id,
                            bnd_position_1 as u32,
                            bnd_strand_1.clone(),
                            bnd_operation_1.clone(),
                            bnd_chromosome_2_id,
                            bnd_position_2 as u32,
                            bnd_strand_2.clone(),
                            bnd_operation_2.clone()
                        )
                    } else {
                        (
                            bnd_chromosome_2_id,
                            bnd_position_2 as u32,
                            bnd_strand_2.clone(),
                            bnd_operation_2.clone(),
                            bnd_chromosome_1_id,
                            bnd_position_1 as u32,
                            bnd_strand_1.clone(),
                            bnd_operation_1.clone()
                        )
                    };
                let graph_operation = SequenceOperation::new(
                    chromosome_1_id,
                    position_1,
                    strand_1,
                    orientation_1,
                    chromosome_2_id,
                    position_2,
                    strand_2,
                    orientation_2,
                    insertion,
                    SequenceOperationVariantTypes::Insertion.clone()
                );
                let variant_record = VariantRecord::new(self.read_id, graph_operation);
                variant_records.push(variant_record);
            }
            prev_alignment_record = curr_alignment_record;
        }
        variant_records
    }

    pub fn is_spliced(&self) -> bool {
        for alignment_record in self.alignment_records.iter() {
            let cs_tag: String;
            if let Some(value) = get_tag_value(&alignment_record.record, "cs") {
                cs_tag = value.to_string();
            } else {
                panic!("Could not find the CS tag.");
            }
            if cs_tag.contains("~") {
                return true;
            }
        }
        false
    }
}
