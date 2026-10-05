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


use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use exacto_core::log_info;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::prelude::*;


pub fn integrate_dna_rna_variants(
    dna_variant_records: &Vec<DNAVariantRecord>,
    rna_variant_records: &Vec<AssembledTranscriptVariantRecord>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    max_exon_offset: u16,
    max_transcript_boundary_offset: u32,
    max_intergenic_distance: u32,
    num_threads: usize
) -> Result<Vec<RNAVariantIntegration>, IntegratorError> {
    // Step 1. Create an index of the DNA variants, whose IDs must be unique
    let mut origins: HashMap<u32, &str> = HashMap::new();
    for dna_variant_record in dna_variant_records.iter() {
        if let Some(origin) = origins.insert(dna_variant_record.variant_id, &dna_variant_record.origin) {
            return Err(IntegratorError::DuplicateDNAVariantId {
                variant_id: dna_variant_record.variant_id,
                origin_1: origin.into(),
                origin_2: dna_variant_record.origin.clone()
            });
        }
    }
    let dna_variant_records_index: DNAVariantIndex = DNAVariantIndex::new(dna_variant_records);

    // Step 2. Sort the exons of every listed reference transcript once
    let mut reference_transcripts: HashMap<&str, TranscriptExons> = HashMap::new();
    for rna_variant_record in rna_variant_records.iter() {
        for reference_transcript_id in split_list(&rna_variant_record.reference_transcript_id) {
            if reference_transcripts.contains_key(reference_transcript_id) {
                continue;
            }
            let reference_transcript: &Transcript = gene_annotator
                .get_transcript(reference_transcript_id)
                .ok_or_else(|| IntegratorError::UnknownTranscript {
                    transcript_id: reference_transcript_id.into(),
                    rna_variant_id: rna_variant_record.variant_id,
                    assembled_transcript_name: rna_variant_record.assembled_transcript_name.clone()
                })?;
            reference_transcripts.insert(reference_transcript_id, TranscriptExons::new(reference_transcript));
        }
    }

    // Step 3. Iterate through each RNA variant and integrate it with any DNA variant(s)
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|error| IntegratorError::ThreadPool { num_threads, reason: error.to_string().into() })?;
    let integrations: Vec<RNAVariantIntegration> = thread_pool.install(|| {
        rna_variant_records
            .par_iter()
            .map(|rna_variant_record| {
                let reference_gene_names: Vec<Box<str>> = split_list(&rna_variant_record.reference_gene_name)
                    .map(Box::from)
                    .collect();
                let reference_transcript_ids: Vec<Box<str>> = split_list(&rna_variant_record.reference_transcript_id)
                    .map(Box::from)
                    .collect();
                let listed_transcripts: Vec<&TranscriptExons> = reference_transcript_ids
                    .iter()
                    .map(|reference_transcript_id| &reference_transcripts[&**reference_transcript_id])
                    .filter(|reference_transcript| !reference_transcript.exons.is_empty())
                    .collect();
                let mut integrated_variant: RNAVariantIntegration = RNAVariantIntegration::new(
                    rna_variant_record.assembled_transcript_name.clone(),
                    reference_gene_names,
                    reference_transcript_ids,
                    rna_variant_record.variant_id
                );
                let rna_positions = [
                    (&*rna_variant_record.chromosome_1, rna_variant_record.position_1, VariantPosition::Position1),
                    (&*rna_variant_record.chromosome_2, rna_variant_record.position_2, VariantPosition::Position2)
                ];

                // Candidates: the DNA variants in or near each listed transcript, and those near
                // each RNA position that lies outside the listed transcripts
                let mut dna_candidates: Vec<&DNAVariantRecord> = Vec::new();
                for listed_transcript in listed_transcripts.iter() {
                    dna_candidates.extend(dna_variant_records_index.get_by_range(
                        listed_transcript.chromosome,
                        listed_transcript.start().saturating_sub(max_transcript_boundary_offset),
                        listed_transcript.end().saturating_add(max_transcript_boundary_offset)
                    ));
                }
                for (chromosome, position, _) in rna_positions.iter() {
                    if is_outside_transcripts(&listed_transcripts, chromosome, *position) {
                        dna_candidates.extend(dna_variant_records_index.get_by_range(
                            chromosome,
                            position.saturating_sub(max_intergenic_distance),
                            position.saturating_add(max_intergenic_distance)
                        ));
                    }
                }
                dna_candidates.sort_unstable_by_key(|dna_variant_record| dna_variant_record.variant_id);
                dna_candidates.dedup_by_key(|dna_variant_record| dna_variant_record.variant_id);

                for dna_variant_record in dna_candidates {
                    let dna_positions = [
                        (&*dna_variant_record.chromosome_1, dna_variant_record.position_1, VariantPosition::Position1),
                        (&*dna_variant_record.chromosome_2, dna_variant_record.position_2, VariantPosition::Position2)
                    ];
                    let mut dna_variant_match: Option<DNAVariantMatch> = None;
                    for (rna_chromosome, rna_position, rna_variant_position) in rna_positions.iter() {
                        for (dna_chromosome, dna_position, dna_variant_position) in dna_positions.iter() {
                            if rna_chromosome != dna_chromosome {
                                continue;
                            }
                            let distance: u32 = rna_position.abs_diff(*dna_position);
                            if dna_variant_match.as_ref().is_some_and(|m| m.get_distance() <= distance) {
                                continue;
                            }
                            if are_positions_linked(
                                &listed_transcripts,
                                rna_chromosome,
                                *rna_position,
                                *dna_position,
                                max_exon_offset,
                                max_transcript_boundary_offset,
                                max_intergenic_distance
                            ) {
                                dna_variant_match = Some(DNAVariantMatch::new(
                                    distance,
                                    rna_variant_position.clone(),
                                    dna_variant_position.clone()
                                ));
                            }
                        }
                    }
                    if let Some(dna_variant_match) = dna_variant_match {
                        integrated_variant.add_dna_variant_id(dna_variant_record.variant_id, dna_variant_match);
                    }
                }

                integrated_variant
            })
            .filter(|integrated_variant| !integrated_variant.dna_variant_ids.is_empty())
            .collect()
    });

    log_info!("Completed integrating DNA and RNA variants.");

    Ok(integrations)
}

/// The exons of a reference transcript in ascending coordinate order.
struct TranscriptExons<'a> {
    chromosome: &'a str,
    exons: Vec<&'a Exon>
}

/// Where a position lies relative to a transcript, by exon number.
#[derive(Debug, PartialEq)]
enum Location {
    Exonic(u16),
    Intronic(u16, u16),

    /// Outside the transcript, `distance` bases beyond the exon at the nearer end.
    Intergenic { exon_number: u16, distance: u32 }
}

impl<'a> TranscriptExons<'a> {
    fn new(transcript: &'a Transcript) -> Self {
        let mut exons: Vec<&Exon> = transcript.exons.values().collect();
        exons.sort_by_key(|exon| exon.start);
        Self { chromosome: &transcript.chromosome, exons }
    }

    fn start(&self) -> u32 {
        self.exons[0].start
    }

    fn end(&self) -> u32 {
        self.exons[self.exons.len() - 1].end
    }

    fn locate(&self, position: u32) -> Location {
        let i: usize = self.exons.partition_point(|exon| exon.start <= position);
        if i == 0 {
            return Location::Intergenic { exon_number: self.exons[0].exon_number, distance: self.exons[0].start - position };
        }
        let exon: &Exon = self.exons[i - 1];
        if position <= exon.end {
            Location::Exonic(exon.exon_number)
        } else if i == self.exons.len() {
            Location::Intergenic { exon_number: exon.exon_number, distance: position - exon.end }
        } else {
            Location::Intronic(exon.exon_number, self.exons[i].exon_number)
        }
    }
}

impl Location {
    /// The exons a position is judged by: its own, the two beside an intron, or the exon at the
    /// end it is near. None when it is farther than `max_transcript_boundary_offset` outside.
    fn exon_numbers(&self, max_transcript_boundary_offset: u32) -> Option<(u16, u16)> {
        match *self {
            Location::Exonic(exon_number) => Some((exon_number, exon_number)),
            Location::Intronic(exon_number_1, exon_number_2) => Some((exon_number_1, exon_number_2)),
            Location::Intergenic { exon_number, distance } => {
                (distance <= max_transcript_boundary_offset).then_some((exon_number, exon_number))
            }
        }
    }
}

fn split_list(list: &str) -> impl Iterator<Item = &str> {
    list.split(LIST_SEPARATOR).filter(|item| !item.is_empty())
}

fn is_outside_transcripts(listed_transcripts: &[&TranscriptExons], chromosome: &str, position: u32) -> bool {
    listed_transcripts
        .iter()
        .filter(|listed_transcript| listed_transcript.chromosome == chromosome)
        .all(|listed_transcript| matches!(listed_transcript.locate(position), Location::Intergenic { .. }))
}

fn are_positions_linked(
    listed_transcripts: &[&TranscriptExons],
    chromosome: &str,
    rna_position: u32,
    dna_position: u32,
    max_exon_offset: u16,
    max_transcript_boundary_offset: u32,
    max_intergenic_distance: u32
) -> bool {
    let mut is_outside: bool = true;
    for listed_transcript in listed_transcripts.iter().filter(|listed_transcript| listed_transcript.chromosome == chromosome) {
        let rna_location: Location = listed_transcript.locate(rna_position);
        let dna_location: Location = listed_transcript.locate(dna_position);
        if let (Some((rna_exon_1, rna_exon_2)), Some((dna_exon_1, dna_exon_2))) = (
            rna_location.exon_numbers(max_transcript_boundary_offset),
            dna_location.exon_numbers(max_transcript_boundary_offset)
        ) {
            if [rna_exon_1, rna_exon_2].iter().any(|rna_exon| {
                [dna_exon_1, dna_exon_2].iter().any(|dna_exon| rna_exon.abs_diff(*dna_exon) <= max_exon_offset)
            }) {
                return true;
            }
        }
        is_outside &= matches!(rna_location, Location::Intergenic { .. }) && matches!(dna_location, Location::Intergenic { .. });
    }
    is_outside && rna_position.abs_diff(dna_position) <= max_intergenic_distance
}


#[cfg(test)]
#[path = "../tests/pipeline/variant_integration.rs"]
mod tests;