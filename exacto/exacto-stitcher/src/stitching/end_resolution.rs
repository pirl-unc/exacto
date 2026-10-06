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
use std::ops::RangeInclusive;

use crate::common::enums::*;
use crate::polyadenylation::signal::{find_polyadenylation_signals, PolyadenylationSignalMatch};
use crate::polyadenylation::tail::{score_polya, PolyATailScore};
use crate::stitching::end_matching::EndMatch;


#[derive(Clone,Debug)]
pub struct EndResolution {
    /// 5' or 3'.
    pub terminus: TranscriptTerminus,

    pub reference_transcript_match: ReferenceTranscriptMatch,

    /// Reference transcript sequence to stitch (add).
    pub stitch_sequence: Box<str>,

    /// For a 5′ end, this is the first retained read base; `stitched_sequence`
    /// is prepended immediately before it.
    /// For a 3′ end, this is the last retained read base; `stitched_sequence`
    /// is appended immediately after it.
    ///
    /// This position uses original read coordinates, before trimming or stitching.
    pub read_join_position: u32
}


/// Resolves the 5' end against every best-ranked end match and keeps it only when they all
/// give the same stitch. Isoforms that share the read's exons but start differently would
/// each stitch a different start, and the read cannot tell which one it came from.
pub(crate) fn resolve_five_prime_end(
    transcript_model: &TranscriptModel,
    end_matches: &[EndMatch],
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> Option<EndResolution> {
    agree(end_matches.iter().map(|end_match| resolve_five_prime_end_match(
        transcript_model,
        end_match,
        gene_annotator,
        chromosome_names_map,
        fasta_map
    )))
}


fn resolve_five_prime_end_match(
    transcript_model: &TranscriptModel,
    end_match: &EndMatch,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap
) -> Option<EndResolution> {
    let bases: &[AlignmentModelBase] = transcript_model.get_alignment_model().get_bases();

    assert_eq!(end_match.terminus, TranscriptTerminus::FivePrime);

    // Step 1. Get the reference transcript sequence.
    let reference_transcript_sequence: ReferenceTranscriptSequence =
        reference_transcript_sequence(
            end_match,
            chromosome_names_map,
            gene_annotator,
            fasta_map
        );
    let reference_bases: &Vec<ReferenceBase> = reference_transcript_sequence.get_bases();

    // Step 2. Get the terminal aligned base.
    let terminal_base: &AlignmentModelBase = &bases[end_match.read_position as usize];

    // Step 3. Contextualize the terminal aligned base in the reference transcript.
    let (chromosome_id, reference_position, strand) = terminal_base
        .get_placement()
        .get_coordinate()
        .unwrap();
    let context: TranscriptTerminusContext = if reference_bases
        .iter()
        .any(|base| base.reference_position == reference_position) {
        TranscriptTerminusContext::Exonic
    } else if reference_bases.windows(2).any(|pair| {
        let low: u32 = pair[0]
            .reference_position
            .min(pair[1].reference_position);
        let high: u32 = pair[0]
            .reference_position
            .max(pair[1].reference_position);
        low < reference_position && reference_position < high
    }) {
        TranscriptTerminusContext::Intronic
    } else {
        TranscriptTerminusContext::Intergenic
    };

    // Step 4. Only an exonic start is stitched. An intronic start could be a novel (cryptic)
    // exon, and an intergenic one lies outside the reference transcript.
    if context != TranscriptTerminusContext::Exonic {
        return None;
    }
    let terminal_reference_base_index: usize = reference_bases
        .iter()
        .position(|base| base.reference_position == reference_position)
        .unwrap();

    // Step 5. Check the read bases upstream of the terminal aligned base (a soft clip).
    let upstream_read_sequence: String = bases[..end_match.read_position as usize]
        .iter()
        .map(|base| base.get_nucleotide().as_str())
        .collect::<String>()
        .to_uppercase();
    let window_start: usize = terminal_reference_base_index.saturating_sub(upstream_read_sequence.len() + 1);
    let upstream_reference_sequence: String = reference_bases[window_start..terminal_reference_base_index]
        .iter()
        .map(|base| base.reference_nucleotide.as_str())
        .collect::<String>()
        .to_uppercase();
    let num_spelled_reference_bases: usize = count_spelled_reference_bases(
        &upstream_read_sequence,
        &upstream_reference_sequence
    )?;

    // Step 6. Stitch the reference bases before them.
    let stitch_sequence: Box<str> = reference_bases[..terminal_reference_base_index - num_spelled_reference_bases]
        .iter()
        .map(|base| base.reference_nucleotide.as_str())
        .collect::<String>()
        .into_boxed_str();
    Some(EndResolution {
        terminus: TranscriptTerminus::FivePrime,
        reference_transcript_match: end_match.reference_transcript_match.clone(),
        stitch_sequence: stitch_sequence,
        read_join_position: 0
    })
}


pub(crate) fn resolve_three_prime_end(
    transcript_model: &TranscriptModel,
    end_matches: &[EndMatch],
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,

    // Polyadenylation hexamer signal
    pas_search_size: usize,
    pas_hexamers: &[&str],
    pas_start_offset_range: &RangeInclusive<usize>,

    // Polyadenylation tail
    polya_tail_min_adenosine_fraction: f64,
    polya_window_size: usize,
    polya_min_consecutive_a: usize
) -> Option<EndResolution> {
    agree(end_matches.iter().map(|end_match| resolve_three_prime_end_match(
        transcript_model,
        end_match,
        gene_annotator,
        chromosome_names_map,
        fasta_map,
        pas_search_size,
        pas_hexamers,
        pas_start_offset_range,
        polya_tail_min_adenosine_fraction,
        polya_window_size,
        polya_min_consecutive_a
    )))
}


fn resolve_three_prime_end_match(
    transcript_model: &TranscriptModel,
    end_match: &EndMatch,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,

    // Polyadenylation hexamer signal
    pas_search_size: usize,
    pas_hexamers: &[&str],
    pas_start_offset_range: &RangeInclusive<usize>,

    // Polyadenylation tail
    polya_tail_min_adenosine_fraction: f64,
    polya_window_size: usize,
    polya_min_consecutive_a: usize
) -> Option<EndResolution> {
    let alignment_model: &AlignmentModel = transcript_model.get_alignment_model();
    let bases: &[AlignmentModelBase] = alignment_model.get_bases();
    let terminal_read_position: usize = end_match.read_position as usize;

    assert_eq!(end_match.terminus, TranscriptTerminus::ThreePrime);

    // Step 1. Get the reference transcript sequence.
    let reference_transcript_sequence: ReferenceTranscriptSequence =
        reference_transcript_sequence(
            end_match,
            chromosome_names_map,
            gene_annotator,
            fasta_map
        );
    let reference_bases: &Vec<ReferenceBase> = reference_transcript_sequence.get_bases();

    // Step 2. Get the terminal aligned base.
    let terminal_base: &AlignmentModelBase = &bases[terminal_read_position];

    // Step 3. Contextualize the terminal aligned base in the reference transcript.
    let (chromosome_id, reference_position, strand) = terminal_base
        .get_placement()
        .get_coordinate()
        .unwrap();
    let chromosome: &Box<str> = chromosome_names_map.get_by_right(&chromosome_id).unwrap();
    let context: TranscriptTerminusContext = if reference_bases
        .iter()
        .any(|base| base.reference_position == reference_position) {
        TranscriptTerminusContext::Exonic
    } else if reference_bases.windows(2).any(|pair| {
        let low: u32 = pair[0]
            .reference_position
            .min(pair[1].reference_position);
        let high: u32 = pair[0]
            .reference_position
            .max(pair[1].reference_position);
        low < reference_position && reference_position < high
    }) {
        TranscriptTerminusContext::Intronic
    } else {
        TranscriptTerminusContext::Intergenic
    };
    if context == TranscriptTerminusContext::Intergenic {
        return None;
    }

    // Step 4. Check the read bases downstream of the terminal aligned base.
    // A polyA tail can be here because the aligner cannot place it, so it arrives as a
    // soft clip, and the join on the terminal aligned base drops it. Anything else
    // (insertions, a non-A soft clip) could be the remnant of a degraded fusion gene
    // or aberrant splicing. In such cases, we do NOT stitch.
    let downstream_read_sequence: String = bases[terminal_read_position + 1..]
        .iter()
        .map(|base| base.get_nucleotide().as_str())
        .collect();
    if !downstream_read_sequence.is_empty()
        && score_polya(&downstream_read_sequence).adenosine_fraction < polya_tail_min_adenosine_fraction {
        return None;
    }

    // Step 5. Look for internal priming: an A-tract in the genome around the terminal
    // aligned base, which an oligo(dT) primer can anneal to.
    // The window spans the base on both sides. This is because:
    // A soft-clipped tail leaves a primed A-tract
    // downstream of it, while a tail that aligned onto the tract leaves the
    // terminal base at the tract's far end, with the As upstream of it.
    // For example:
    // genome      ...G C T G A A A A A A A A A A G T C C T G...
    // read        ...G C T G A A A A A A A A A A
    //                                          ^ terminal aligned base
    //             upstream: the whole tract    | downstream: G T C C T G
    let start: usize = (reference_position as usize)
        .saturating_sub(polya_window_size)
        .max(1);
    let end: usize = reference_position as usize + polya_window_size;
    let flanking_genomic_sequence: String = fasta_map
        .try_get_sequence(chromosome, start, end)
        .unwrap()
        .to_uppercase();

    // The genomic interval is the same on both strands. Only the orientation
    // differs: As on the transcript are Ts on the genome for a reverse-strand
    // gene, so the window is reverse complemented into transcript orientation.
    let flanking_reference_sequence: Box<str> = if *strand == Strand::Forward {
        flanking_genomic_sequence.into_boxed_str()
    } else {
        reverse_complement(&flanking_genomic_sequence)
    };
    let polya_score: PolyATailScore = score_polya(&*flanking_reference_sequence);
    let is_internally_primed: bool = polya_score.max_consecutive_adenosine >= polya_min_consecutive_a;

    // Step 6. Resolve the end.
    match context {
        TranscriptTerminusContext::Exonic => {
            // A polyA tail the genome does not template, after a polyadenylation signal,
            // marks a polyadenylation site: the transcript ends here although the reference
            // transcript goes on (alternative polyadenylation). Therefore, do NOT stitch.
            if !downstream_read_sequence.is_empty()
                && !is_internally_primed
                && has_polyadenylation_signal(bases, terminal_read_position, pas_search_size, pas_hexamers, pas_start_offset_range) {
                return None;
            }

            // Decide which reference bases need to be added.
            let terminal_reference_base_index: usize = reference_bases
                .iter()
                .position(|base| base.reference_position == reference_position)
                .unwrap();
            let stitch_sequence: Box<str> = reference_bases[terminal_reference_base_index + 1..]
                .iter()
                .map(|base| base.reference_nucleotide.as_str())
                .collect::<String>()
                .into_boxed_str();
            Some(EndResolution {
                terminus: TranscriptTerminus::ThreePrime,
                reference_transcript_match: end_match.reference_transcript_match.clone(),
                stitch_sequence: stitch_sequence,
                read_join_position: end_match.read_position
            })
        },
        TranscriptTerminusContext::Intronic => {
            // Only internal priming shows that the end is an artifact. Without it the read
            // could end at an intronic polyadenylation site, or be a nascent or
            // intron-retaining RNA, and rolling back would delete its intronic bases.
            // Therefore, do NOT stitch.
            if !is_internally_primed {
                return None;
            }

            // The end is an artifact.
            // Roll back to the last base of the exon before the intron.
            let target_reference_base_index: usize = reference_bases.windows(2).position(|pair| {
                let low: u32 = pair[0].reference_position.min(pair[1].reference_position);
                let high: u32 = pair[0].reference_position.max(pair[1].reference_position);
                low < end_match.reference_position && end_match.reference_position < high
            })?;
            let target_reference_position: u32 = reference_bases[target_reference_base_index].reference_position;

            // A read can cross the same locus more than once; the crossing nearest the 3' end wins.
            let target_read_position: u32 = bases
                .iter()
                .filter(|base| matches!(
                    base.get_kind(),
                    AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch
                ))
                .filter(|base| {
                    base.get_placement().get_coordinate() == Some((end_match.chromosome_id, target_reference_position, &end_match.strand))
                })
                .map(|base| base.get_read_position())
                .max()?;

            // Refuse a roll-back that trims across an alignment event: the trimmed
            // stretch then holds structure the reference does not know (e.g. a novel exon).
            let crosses_event: bool = alignment_model.get_events().values().any(|event| {
                event.get_prev_read_position() >= target_read_position
                    && event.get_prev_read_position() < end_match.read_position
            });
            if crosses_event {
                return None;
            }

            Some(EndResolution {
                terminus: TranscriptTerminus::ThreePrime,
                reference_transcript_match: end_match.reference_transcript_match.clone(),
                stitch_sequence: reference_bases[target_reference_base_index + 1..]
                    .iter()
                    .map(|base| base.reference_nucleotide.as_str())
                    .collect::<String>()
                    .into_boxed_str(),
                read_join_position: target_read_position
            })
        },
        TranscriptTerminusContext::Intergenic => {
            None
        }
    }
}


/// Whether the read bases ending at the terminal aligned base hold a polyadenylation signal
/// hexamer in range.
fn has_polyadenylation_signal(
    bases: &[AlignmentModelBase],
    terminal_read_position: usize,
    pas_search_size: usize,
    pas_hexamers: &[&str],
    pas_start_offset_range: &RangeInclusive<usize>
) -> bool {
    let end: usize = terminal_read_position + 1;
    let start: usize = end.saturating_sub(pas_search_size);
    let upstream_read_sequence: String = bases[start..end]
        .iter()
        .map(|base| base.get_nucleotide().as_str())
        .collect();
    let pas_matches: Vec<PolyadenylationSignalMatch> = find_polyadenylation_signals(
        &upstream_read_sequence,
        pas_hexamers,
        pas_start_offset_range
    );
    !pas_matches.is_empty()
}


/// The resolution every end match gives, when they all give the same stitch at the same
/// join; `None` when they differ, when any one of them declines, or when there is none.
fn agree(mut resolutions: impl Iterator<Item = Option<EndResolution>>) -> Option<EndResolution> {
    let first: EndResolution = resolutions.next()??;
    for resolution in resolutions {
        let resolution: EndResolution = resolution?;
        if !resolution.stitch_sequence.eq_ignore_ascii_case(&first.stitch_sequence)
            || resolution.read_join_position != first.read_join_position {
            return None;
        }
    }
    Some(first)
}


/// How many reference bases, ending right before the terminal aligned base, the read bases
/// before it spell: of the reference suffixes one base shorter to one base longer than the
/// read bases, the one at the fewest edits, if that is at most one. `None` when no suffix is
/// that close, or when two are equally close, so the boundary is ambiguous.
/// `upstream_reference_sequence` ends right before the terminal aligned base.
fn count_spelled_reference_bases(
    upstream_read_sequence: &str,
    upstream_reference_sequence: &str
) -> Option<usize> {
    let read: &[u8] = upstream_read_sequence.as_bytes();
    let reference: &[u8] = upstream_reference_sequence.as_bytes();
    let mut best: Option<(u8, usize)> = None;
    let mut is_tied: bool = false;
    for length in read.len().saturating_sub(1)..=(read.len() + 1).min(reference.len()) {
        let Some(num_edits) = count_edits_up_to_one(read, &reference[reference.len() - length..]) else {
            continue;
        };
        match best {
            Some((best_num_edits, _)) if num_edits > best_num_edits => {},
            Some((best_num_edits, _)) if num_edits == best_num_edits => is_tied = true,
            _ => {
                best = Some((num_edits, length));
                is_tied = false;
            }
        }
    }
    if is_tied {
        return None;
    }
    best.map(|(_, length)| length)
}


/// 0 when `a` equals `b`, 1 when one substitution, insertion or deletion turns one into the
/// other, `None` otherwise. One edit, if any, sits at the first position where they differ.
fn count_edits_up_to_one(a: &[u8], b: &[u8]) -> Option<u8> {
    let prefix_length: usize = a.iter().zip(b).take_while(|(x, y)| x == y).count();
    let (a, b) = (&a[prefix_length..], &b[prefix_length..]);
    if a.is_empty() && b.is_empty() {
        return Some(0);
    }
    let is_one_edit: bool = a.get(1..) == Some(b)
        || b.get(1..) == Some(a)
        || (a.len() == b.len() && a[1..] == b[1..]);
    is_one_edit.then_some(1)
}


fn reference_transcript_sequence(
    end_match: &EndMatch,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    gene_annotator: &impl GeneAnnotator,
    fasta_map: &FastaMap
) -> ReferenceTranscriptSequence {
    let transcript_id: &str = end_match.reference_transcript_match.get_reference_transcript_id();
    let transcript: &Transcript = gene_annotator
        .get_transcript(transcript_id)
        .unwrap_or_else(|| panic!("End match names transcript {transcript_id}, which the annotator does not have."));
    ReferenceTranscriptSequence::from_reference_transcript(transcript, chromosome_names_map, fasta_map)
}


#[cfg(test)]
#[path = "../tests/stitching/end_resolution.rs"]
mod tests;