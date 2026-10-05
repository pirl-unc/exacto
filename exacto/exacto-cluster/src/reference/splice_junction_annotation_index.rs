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
use exacto_caller::prelude::SpliceJunction;
use exacto_core::prelude::{GeneAnnotator, Intron, Strand};
use std::collections::{BTreeSet, HashMap, HashSet};


pub struct SpliceJunctionAnnotationIndex {
    /// HashMap<(chromosome ID, intron start, intron end), BTreeSet<gene ID>>
    pub junctions_to_genes: HashMap<(u16, u32, u32), BTreeSet<Box<str>>>,

    /// The pairs of introns that follow one another in a reference transcript, in transcript
    /// order.
    pub annotated_transitions: HashSet<(SpliceJunction, SpliceJunction)>,

    /// The full span of every annotated gene.
    /// HashMap<chromosome ID, HashMap<gene ID, (start, end)>>
    pub gene_spans: HashMap<u16, HashMap<Box<str>, (u32, u32)>>,

    /// (gene ID, transcript ID) of every reference transcript that holds an intron.
    pub transcripts: Vec<(Box<str>, Box<str>)>,

    /// HashMap<(chromosome ID, intron start, intron end, strand), Vec<index into `transcripts`>>
    pub junctions_to_transcripts: HashMap<(u16, u32, u32, Strand), Vec<u32>>
}

impl SpliceJunctionAnnotationIndex {
    pub fn new(
        gene_annotator: &impl GeneAnnotator,
        chromosome_names_map: &BiMap<Box<str>, u16>
    ) -> Self {
        let mut junctions_to_genes: HashMap<(u16, u32, u32), BTreeSet<Box<str>>> = HashMap::new();
        let mut annotated_transitions: HashSet<(SpliceJunction, SpliceJunction)> = HashSet::new();
        let mut gene_spans: HashMap<u16, HashMap<Box<str>, (u32, u32)>> = HashMap::new();
        let mut transcripts: Vec<(Box<str>, Box<str>)> = Vec::new();
        let mut junctions_to_transcripts: HashMap<(u16, u32, u32, Strand), Vec<u32>> = HashMap::new();
        for transcript in gene_annotator.get_transcripts() {
            // Step 1. The span of the gene.
            if let Some(&chromosome_id) = chromosome_names_map.get_by_left(&transcript.chromosome) {
                let start: u32 = transcript.start.min(transcript.end);
                let end: u32 = transcript.start.max(transcript.end);
                gene_spans
                    .entry(chromosome_id)
                    .or_default()
                    .entry(transcript.gene_id.clone())
                    .and_modify(|span| {
                        span.0 = span.0.min(start);
                        span.1 = span.1.max(end);
                    })
                    .or_insert((start, end));
            }

            // Step 2. The introns, and the pairs of them in transcript order.
            let introns: Vec<Intron> = transcript.get_introns();
            let mut reference_chain: Vec<SpliceJunction> = Vec::new();
            for intron in introns.iter() {
                let Some(&chromosome_id) = chromosome_names_map.get_by_left(&intron.chromosome) else {
                    continue;
                };
                let start: u32 = intron.start.min(intron.end);
                let end: u32 = intron.start.max(intron.end);
                junctions_to_genes
                    .entry((chromosome_id, start, end))
                    .or_default()
                    .insert(intron.gene_id.clone());
                junctions_to_transcripts
                    .entry((chromosome_id, start, end, intron.strand.clone()))
                    .or_default()
                    .push(transcripts.len() as u32);
                let (donor, acceptor): (u32, u32) = if intron.strand == Strand::Reverse {
                    (intron.end, intron.start)
                } else {
                    (intron.start, intron.end)
                };
                reference_chain.push(SpliceJunction::new(
                    chromosome_id,
                    chromosome_id,
                    donor,
                    acceptor,
                    intron.strand.clone(),
                    intron.strand.clone()
                ));
            }
            for pair in reference_chain.windows(2) {
                annotated_transitions.insert((pair[0].clone(), pair[1].clone()));
            }
            if let Some(intron) = introns.first() {
                transcripts.push((intron.gene_id.clone(), intron.transcript_id.clone()));
            }
        }

        Self {
            junctions_to_genes,
            annotated_transitions,
            gene_spans,
            transcripts,
            junctions_to_transcripts
        }
    }
}


#[cfg(test)]
#[path = "../tests/reference/splice_junction_annotation_index.rs"]
mod tests;