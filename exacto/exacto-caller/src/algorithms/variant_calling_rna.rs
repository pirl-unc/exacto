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
use exacto_util::structs::gene_annotation::transcript::Transcript;
use exacto_util::prelude::*;
use noodles_bam as bam;
use rayon::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use std::collections::{HashMap, HashSet};

use crate::prelude::*;
use crate::structs::transcript_model_set::TranscriptModelSet;


pub fn identify_variant_transcripts(
    bam_file: &str,
    bam_bai_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    min_mapping_quality: usize,
    min_average_base_quality: f32,
    num_threads: usize
) -> TranscriptModelSet {
    // Step 1. Get a map of read names and IDs
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        num_threads
    );

    // Step 2. Get chromosome names map
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);

    // Step 3. Fetch all BAM records
    let records_map: HashMap<usize,Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        num_threads
    );

    // Step 4. Construct transcript models
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let variant_transcript_models: Vec<TranscriptModel> = thread_pool.install(|| {
        records_map
            .par_iter()
            .filter_map(|(read_id, records)| {
                let read_sequence: Box<str> =  get_original_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
                let quality_scores: Vec<u8> = get_original_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
                let alignment: Alignment = Alignment::new(
                    *read_id,
                    read_sequence,
                    quality_scores,
                    records.clone()
                );

                // Identify exons
                let exons: Vec<TranscriptModelExon> = alignment.identify_exons(min_mapping_quality);

                // Identify splice junctions
                let splice_junctions: Vec<TranscriptModelSpliceJunction> = alignment.identify_splice_junctions(
                    &chromosome_names_map,
                    reference_genome_fasta_file,
                    min_mapping_quality
                );

                // Identify overlapping gene IDs
                let mut reference_gene_ids: HashSet<Box<str>> = Alignment::identify_overlapping_gene_ids(
                    &exons,
                    gene_annotator,
                    &chromosome_names_map
                );

                // Identify closest reference transcript IDs
                let reference_transcript_ids: Vec<Box<str>> = Alignment::identify_closest_reference_transcript_ids(
                    &exons,
                    gene_annotator,
                    &chromosome_names_map,
                    &reference_gene_ids
                );

                // Identify splice variant records
                let splice_variant_records: Vec<VariantRecord> = alignment.identify_splice_variant_records(
                    gene_annotator,
                    &chromosome_names_map,
                    min_mapping_quality
                );

                // Identify sequence variant records
                let sequence_variant_records: Vec<VariantRecord> = alignment.identify_sequence_variant_records(
                    min_mapping_quality,
                    min_average_base_quality
                );

                // Construct a transcript model
                let read_ids: Vec<usize> = vec![*read_id];
                let mut transcript_model: TranscriptModel = TranscriptModel::new(
                    alignment.read_id,
                    reference_transcript_ids,
                    read_ids
                );
                for exon in exons {
                    transcript_model.add_exon(exon);
                }
                for splice_junction in splice_junctions {
                    transcript_model.add_splice_junction(splice_junction);
                }
                for splice_variant_record in splice_variant_records {
                    let mut variant_call: VariantCall = VariantCall::new();
                    variant_call.add_variant_record(splice_variant_record);
                    transcript_model.add_variant_call(variant_call);
                }
                for sequence_variant_record in sequence_variant_records {
                    let mut variant_call: VariantCall = VariantCall::new();
                    variant_call.add_variant_record(sequence_variant_record);
                    transcript_model.add_variant_call(variant_call);
                }

                if is_nascent_transcript(&transcript_model, gene_annotator, &chromosome_names_map) {
                    None
                } else {
                    if transcript_model.variant_calls.is_empty() == false {
                        Some(transcript_model)
                    } else {
                        None
                    }
                }
            })
            .collect()
    });

    let mut transcript_model_set: TranscriptModelSet = TranscriptModelSet::new();
    let mut transcript_id: usize = 1;;
    for mut variant_transcript_model in variant_transcript_models {
        variant_transcript_model.transcript_id = transcript_id;
        transcript_model_set.add_transcript_model(variant_transcript_model);
        transcript_id += 1;
    }

    transcript_model_set.load_read_names(read_names_map);
    transcript_model_set.load_chromosome_names(chromosome_names_map);

    transcript_model_set
}

pub fn is_nascent_transcript(
    transcript_model: &TranscriptModel,
    gene_annotator: &impl GeneAnnotator,
    chromosome_names_map: &BiMap<Box<str>,u16>
) -> bool {
    if transcript_model.splice_junctions.is_empty() {
        // Identify all relevant gene IDs
        let mut gene_ids: HashSet<Box<str>> = HashSet::new();
        for exon in transcript_model.exons.iter() {
            let chromosome: Box<str> = chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap().to_string().into_boxed_str();
            let gene_ids_: Vec<Box<str>> = gene_annotator.get_gene_ids_overlapping_region(&*chromosome, exon.start, exon.end);
            for gene_id in gene_ids_ {
                gene_ids.insert(gene_id);
            }
        }

        // This transcript is not a nascent transcript if it overlaps with a transcript with only 1 exon
        for gene_id in gene_ids.iter() {
            let gene: &Gene = gene_annotator.get_gene(gene_id).unwrap();
            for transcript_id in gene.get_transcript_ids().iter() {
                let transcript: &Transcript = gene_annotator.get_transcript(transcript_id).unwrap();
                if transcript.get_exon_ids().len() == 1 {
                    let exon_: &Exon = transcript.get_exon(transcript.get_exon_ids().first().unwrap()).unwrap();
                    for exon in transcript_model.exons.iter() {
                        if overlaps(exon.start as isize, exon.end as isize, exon_.start as isize, exon_.end as isize) {
                            return false;
                        }
                    }
                }
            }
        }
        true
    } else {
        false
    }
}
