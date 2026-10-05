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
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use std::cmp::Reverse;


pub(crate) struct EndMatch {
    pub terminus: TranscriptTerminus,
    pub reference_transcript_match: ReferenceTranscriptMatch,

    /// Based on AlignmentModelBase
    pub chromosome_id: u16,
    pub reference_position: u32,
    pub strand: Strand,
    pub read_position: u32
}


pub(crate) fn identify_transcript_terminus_matches(
    transcript_model: &TranscriptModel,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    terminus: &TranscriptTerminus,
    min_num_splice_junction_matches: usize
) -> Vec<EndMatch> {
    // Step 1. Locate the terminal aligned base.
    let Some((chromosome_id, reference_position, strand, read_position)) = locate_terminal_aligned_base(
        transcript_model.get_alignment_model(),
        terminus
    ) else {
        return Vec::new();
    };

    let Some(chromosome) = chromosome_names_map.get_by_right(&chromosome_id) else {
        return Vec::new();
    };

    // Step 2. Score every transcript of every gene the read matched.
    let mut gene_ids: Vec<&str> = transcript_model
        .get_reference_transcript_matches()
        .iter()
        .map(|m| m.get_reference_gene_id())
        .collect();
    gene_ids.sort();
    gene_ids.dedup();

    let mut candidates: Vec<ReferenceTranscriptMatch> = Vec::new();
    for gene in gene_ids.iter().filter_map(|gene_id| gene_annotator.get_gene(gene_id)) {
        for transcript in gene.transcripts.values() {
            if transcript.chromosome != *chromosome
                || transcript.strand != strand
                || reference_position < transcript.start
                || transcript.end < reference_position {
                continue;
            }
            let m: ReferenceTranscriptMatch = score_reference_transcript(
                transcript_model.get_exons(),
                transcript_model.get_splice_junctions(),
                transcript,
                gene,
                chromosome_names_map
            );

            // A single-exon reference transcript has no splice junctions to
            // match, so it bypasses the splice-junction minimum.
            if transcript.exons.len() == 1
                || m.num_splice_junction_matches() >= min_num_splice_junction_matches {
                candidates.push(m);
            }
        }
    }

    // Step 3. Rank candidates. Lower tuples win:
    // 1. More matched splice junctions.
    // 2. Fewer terminal query-only bases.
    // 3. Fewer internal query-only bases.
    // 4. More overlapping bases.
    // Keep every candidate with the best score.
    let score = |m: &ReferenceTranscriptMatch| {
        (
            Reverse(m.num_splice_junction_matches()),
            m.num_terminal_query_only_bases(),
            m.num_internal_query_only_bases(),
            Reverse(m.num_overlapping_bases())
        )
    };
    let Some(best_score) = candidates.iter().map(|m| score(m)).min() else {
        return Vec::new();
    };
    candidates.retain(|m| score(m) == best_score);

    // Step 4. Order the best candidates by the annotator's rank, from transcript ID order, so
    // that the reported transcript never depends on hash order.
    candidates.sort_by(|a, b| a.get_reference_transcript_id().cmp(b.get_reference_transcript_id()));
    let transcripts: Vec<&Transcript> = gene_annotator.rank_transcripts(
        candidates
            .iter()
            .filter_map(|m| gene_annotator.get_transcript(m.get_reference_transcript_id()))
            .collect()
    );

    // Step 5. GENCODE tags a transcript whose 5' or 3' end was not found. Its annotated end
    // there is where the evidence stopped, so it has no end to copy: decline.
    let incomplete_end_tag: &str = match terminus {
        TranscriptTerminus::FivePrime => "mRNA_start_NF",
        TranscriptTerminus::ThreePrime => "mRNA_end_NF"
    };
    if transcripts.iter().any(|transcript| transcript.tags.contains(incomplete_end_tag)) {
        return Vec::new();
    }

    transcripts
        .iter()
        .filter_map(|transcript| candidates
            .iter()
            .find(|m| m.get_reference_transcript_id() == &*transcript.transcript_id))
        .map(|m| EndMatch {
            terminus: *terminus,
            reference_transcript_match: m.clone(),
            chromosome_id,
            reference_position,
            strand: strand.clone(),
            read_position
        })
        .collect()
}


fn locate_terminal_aligned_base(
    alignment_model: &AlignmentModel,
    terminus: &TranscriptTerminus
) -> Option<(u16, u32, Strand, u32)> {
    let is_anchoring = |base: &&AlignmentModelBase| matches!(
        base.get_kind(),
        AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch
    );

    let base: &AlignmentModelBase = match terminus {
        TranscriptTerminus::FivePrime => alignment_model.get_bases().iter().find(is_anchoring),
        TranscriptTerminus::ThreePrime => alignment_model.get_bases().iter().rev().find(is_anchoring)
    }?;

    let (chromosome_id, position, strand) = base
        .get_placement()
        .get_coordinate()?;

    Some((chromosome_id, position, strand.clone(), base.get_read_position()))
}


#[cfg(test)]
#[path = "../tests/stitching/end_matching.rs"]
mod tests;