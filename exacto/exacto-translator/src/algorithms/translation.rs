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


use exacto_caller::prelude::{DNAVariantRecord, RNAVariantRecord, TranscriptModelStructureRecord};
use exacto_integrator::prelude::IntegratedVariantRecord;
use std::collections::HashSet;

use crate::prelude::*;


pub fn translate_structures(
    assembled_transcript_support_records: &Vec<AssembledTranscriptSupportRecord>,
    transcript_model_structure_records: &Vec<TranscriptModelStructureRecord>,
    rna_variant_records: &Vec<RNAVariantRecord>,
    dna_variant_records: &Vec<DNAVariantRecord>,
    integrated_variant_records: &Vec<IntegratedVariantRecord>,
    translation_strategy: TranslationStrategy,
    start_codons: HashSet<&str>,
    num_threads: usize
) -> TranscriptSet {
    let mut ts: TranscriptSet = build_transcript_set(
        assembled_transcript_support_records,
        transcript_model_structure_records,
        rna_variant_records,
        dna_variant_records,
        integrated_variant_records
    );
    
    ts.translate(
        translation_strategy,
        start_codons,
        num_threads
    );
    
    ts
}

pub fn translate_sequences(
        sequences: Vec<(Box<str>, Box<str>)>,
        translation_strategy: TranslationStrategy,
        start_codons: HashSet<&str>,
        num_threads: usize
) -> TranscriptSet {
    let mut transcripts: Vec<Transcript> = Vec::with_capacity(sequences.len());
    for (transcript_id, sequence) in sequences.into_iter() {
        let transcript: Transcript = Transcript::new_with_default(transcript_id, sequence);
        transcripts.push(transcript);
    }

    let mut ts: TranscriptSet = TranscriptSet::new(transcripts);

    ts.translate(
        translation_strategy,
        start_codons,
        num_threads
    );

    ts
}
