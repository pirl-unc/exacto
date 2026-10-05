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
use std::collections::HashMap;
use std::sync::Arc;

use crate::prelude::*;


pub fn characterize_template_switch(
    variant_call: &VariantCall,
    transcript_models_map: &HashMap<ReadID, Arc<TranscriptModel>>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    flank: u32,
    max_loop_length: u32,
    slack: u32
) -> TemplateSwitchEvidence {
    let graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
    let chromosome_1 = chromosome_names_map.get_by_right(&graph_operation.get_chromosome_1()).unwrap();
    let chromosome_2 = chromosome_names_map.get_by_right(&graph_operation.get_chromosome_2()).unwrap();

    // Step 1. Panel A: flank homology across the junction.
    let junction_homology: Option<JunctionHomology> = compute_junction_homology(
        graph_operation,
        chromosome_names_map,
        fasta_map,
        flank
    );

    // Step 2. Panel B: the inverted repeat at a fold-back. `None` for any other shape.
    let foldback_stem: Option<FoldbackStem> = compute_foldback_stem(
        graph_operation,
        chromosome_names_map,
        fasta_map,
        flank,
        max_loop_length
    );

    // Step 3. Both breakpoints on an observed exon boundary.
    fn is_annotated_exon_boundary(
        gene_annotator: &impl GeneAnnotator,
        chromosome: &str,
        position: ReferencePosition
    ) -> bool {
        gene_annotator
            .get_transcript_ids_overlapping_region(chromosome, position, position)
            .iter()
            .filter_map(|transcript_id| gene_annotator.get_transcript(transcript_id))
            .flat_map(|transcript| transcript.exons.values())
            .any(|exon| exon.start == position || exon.end == position)
    }
    let at_exon_boundaries: bool = is_junction_variant_type(graph_operation.get_variant_type())
        && is_annotated_exon_boundary(gene_annotator, chromosome_1, graph_operation.get_position_1())
        && is_annotated_exon_boundary(gene_annotator, chromosome_2, graph_operation.get_position_2());

    // Step 4. Fold-backs only: does any supporting read leave the mirrored footprint?
    let has_exit_junction: bool = foldback_stem.is_some()
        && variant_call
        .get_read_ids()
        .iter()
        .any(|read_id| read_leaves_footprint(&transcript_models_map[read_id], graph_operation, slack));

    // Step 5. Pooled breakpoints only: member spread beyond the homology interval, inside
    // which a breakpoint is genuinely ambiguous.
    let homology_total: u32 = junction_homology.map_or(0, |h| h.total());
    let dispersion_beyond_homology: Option<u32> = pooled_breakpoint_spread(variant_call)
        .map(|spread| spread.saturating_sub(homology_total));

    TemplateSwitchEvidence {
        junction_homology,
        foldback_stem,
        at_exon_boundaries,
        has_exit_junction,
        dispersion_beyond_homology,
        num_reads: variant_call.get_read_ids().len() as ReadSupport
    }
}


/// The pure decision. Panel B is judged first, on geometry; panel A after, on homology.
/// Every rule of the filter lives here and nowhere else.
pub fn classify_template_switch(
    evidence: &TemplateSwitchEvidence,
    min_homology: u32,
    soft_min_homology: u32,
    max_breakpoint_dispersion: u32,
    foldback_max_distance: u32,
    foldback_min_stem: u32
) -> TemplateSwitchVerdict {
    use TemplateSwitchReason::*;
    use TemplateSwitchVerdict::*;

    // Panel B. A fold-back is decided on its shape and its stem, never on h: the two
    // flanks are one locus read both ways, so h only measures the stem against itself.
    match classify_foldback(evidence, foldback_max_distance, foldback_min_stem) {
        FoldbackVerdict::Hairpin => return Flagged(Foldback),
        FoldbackVerdict::Duplication | FoldbackVerdict::Unresolved => return Clear,
        FoldbackVerdict::NotFoldback => {}
    }

    // Panel A. Two templates sharing a motif at the junction.
    let Some(homology) = evidence.junction_homology else {
        return NotAssessed;
    };
    // The spliceosome could have made this junction, or a rearrangement of exons did.
    // Either way the repeat is not RT's: RT never saw the intron the flank runs into.
    if homology.canonical_splice || evidence.at_exon_boundaries {
        return Clear;
    }
    let h: u32 = homology.total();
    if h >= min_homology {
        return Flagged(StrongHomology);
    }
    if h >= soft_min_homology {
        if evidence.dispersion_beyond_homology.is_some_and(|d| d > max_breakpoint_dispersion) {
            return Flagged(SoftHomologyWithDispersion);
        }
    }
    Clear
}


#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum FoldbackVerdict {
    /// Not fold-back shaped, or farther apart than `foldback_max_distance`.
    NotFoldback,
    
    /// Fold on an exon boundary, or a supporting read leaves the mirrored arm: the shape
    /// of an inverted duplication, not a hairpin.
    Duplication,
    
    /// An inverted repeat long enough to prime sits at the fold: an RT hairpin.
    Hairpin,
    
    /// Fold-back shaped with neither signature: report, do not suppress.
    Unresolved
}


#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct FoldbackStem {
    /// Bases between the two breakpoints: the hairpin's footprint (stem length + loop length + stem length).
    pub span: u32,

    /// Longest run of complementary pairs closing across the fold.
    pub stem: u32,

    /// Unpaired bases between the two halves of that stem.
    pub loop_length: u32
}


#[derive(Clone,Copy,Debug,Eq,PartialEq,Hash)]
pub(crate) struct JunctionHomology {
    pub h_5prime: u32,
    pub h_3prime: u32,
    pub canonical_splice: bool
}

impl JunctionHomology {
    pub fn total(&self) -> u32 {
        self.h_5prime + self.h_3prime
    }
}


struct JunctionSide<'a> {
    pub chromosome: &'a str,
    pub position: ReferencePosition,
    pub strand: Strand,
    pub operation: GraphOperationType
}


enum SideRole {
    Leaves,
    Enters
}


fn classify_foldback(
    evidence: &TemplateSwitchEvidence,
    foldback_max_distance: u32,
    min_stem: u32
) -> FoldbackVerdict {
    let Some(stem) = evidence.foldback_stem else {
        return FoldbackVerdict::NotFoldback;
    };
    
    if stem.span > foldback_max_distance {
        return FoldbackVerdict::NotFoldback;
    }
    
    // A hairpin folds where a stem is, never on an exon end, and never continues into
    // sequence the first arm did not cover: its read ends inside the mirrored arm.
    if evidence.at_exon_boundaries || evidence.has_exit_junction {
        return FoldbackVerdict::Duplication;
    }
    
    if stem.stem >= min_stem {
        return FoldbackVerdict::Hairpin;
    }
    
    FoldbackVerdict::Unresolved
}


/// Classify whether the read leaves or enters the template at this side:
/// On the `Forward` strand, `Downstream` leaves.
/// On the `Forward` strand, `Upstream` enters.
/// On the `Reverse` strand, `Downstream` enters.
/// On the `Reverse` strand, `Upstream` leaves.
fn classify_side_role(side: &JunctionSide) -> Option<SideRole> {
    match (&side.operation, &side.strand) {
        (GraphOperationType::Downstream, Strand::Forward) => Some(SideRole::Leaves),
        (GraphOperationType::Upstream, Strand::Forward) => Some(SideRole::Enters),
        (GraphOperationType::Downstream, Strand::Reverse) => Some(SideRole::Enters),
        (GraphOperationType::Upstream, Strand::Reverse) => Some(SideRole::Leaves),
        _ => None
    }
}


/// Count how many leading bytes two strings share, stopping at the first mismatch or the first N.
fn common_prefix(a: &str, b: &str) -> u32 {
    a.bytes()
        .zip(b.bytes())
        .take_while(|(x, y)| x == y && *x != b'N')
        .count() as u32
}


/// Count how many trailing bytes two strings share, stopping at the first mismatch or the first N
/// when walking backward from the end.
fn common_suffix(a: &str, b: &str) -> u32 {
    a.bytes()
        .rev()
        .zip(b.bytes().rev())
        .take_while(|(x, y)| x == y && *x != b'N')
        .count() as u32
}


fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other   // N never pairs
    }
}


fn compute_foldback_stem(
    graph_operation: &GraphOperation,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    flank: u32,
    max_loop_length: u32
) -> Option<FoldbackStem> {
    if !is_foldback_shape(graph_operation) {
        return None;
    }

    let chromosome: &str = chromosome_names_map.get_by_right(&graph_operation.get_chromosome_1())?;
    let low: ReferencePosition = graph_operation.get_position_1().min(graph_operation.get_position_2());
    let high: ReferencePosition = graph_operation.get_position_1().max(graph_operation.get_position_2());
    let span: u32 = high - low;

    // Reference window: the footprint plus `flank` bases either side, on the + strand.
    // Which strand the read used does not matter: an inverted repeat is one on both.
    let start: ReferencePosition = low.saturating_sub(flank).max(1);
    let end: ReferencePosition = high + flank;
    let window: String = fasta_map
        .try_get_sequence(chromosome, start as usize, end as usize)?
        .to_uppercase();
    let bytes: &[u8] = window.as_bytes();
    let lo: usize = (low - start) as usize;
    let hi: usize = (high - start) as usize;

    let mut best: (u32, u32) = (0, 0);   // (stem, loop)
    for centre in lo..=hi + 1 {
        for loop_length in 0..=(max_loop_length as usize) {
            if centre + loop_length > bytes.len() {
                break;
            }
            let mut stem: usize = 0;
            while centre >= stem + 1 && centre + loop_length + stem < bytes.len() {
                let left: u8 = bytes[centre - 1 - stem];
                let right: u8 = bytes[centre + loop_length + stem];
                if left == b'N' || complement(left) != right {
                    break;
                }
                stem += 1;
            }
            if stem as u32 > best.0 {
                best = (stem as u32, loop_length as u32);
            }
        }
    }

    Some(FoldbackStem { span, stem: best.0, loop_length: best.1 })
}


fn compute_junction_homology(
    graph_operation: &GraphOperation,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    flank: u32
) -> Option<JunctionHomology> {
    if !is_junction_variant_type(graph_operation.get_variant_type()) {
        return None;
    }

    if flank == 0 {
        return None;
    }

    // Step 1. Get each junction side.
    let chromosome_1: &str = chromosome_names_map.get_by_right(&graph_operation.get_chromosome_1()).unwrap();
    let chromosome_2: &str = chromosome_names_map.get_by_right(&graph_operation.get_chromosome_2()).unwrap();
    let side_1: JunctionSide = JunctionSide {
        chromosome: chromosome_1,
        position: graph_operation.get_position_1(),
        strand: graph_operation.get_strand_1().clone(),
        operation: graph_operation.get_operation_type_1().clone()
    };
    let side_2: JunctionSide = JunctionSide {
        chromosome: chromosome_2,
        position: graph_operation.get_position_2(),
        strand: graph_operation.get_strand_2().clone(),
        operation: graph_operation.get_operation_type_2().clone()
    };

    // Step 2. Classify the role of each junction side.
    let role_1: SideRole = classify_side_role(&side_1)?;
    let role_2: SideRole = classify_side_role(&side_2)?;
    let (leaves, enters) = match (role_1, role_2) {
        (SideRole::Leaves, SideRole::Enters) => (side_1, side_2),
        (SideRole::Enters, SideRole::Leaves) => (side_2, side_1),
        _ => return None
    };

    // Step 3. Get the side flank sequences.
    let (a_in, a_out) = get_side_flank_sequences(&leaves, SideRole::Leaves, fasta_map, flank)?;
    let (b_in, b_out) = get_side_flank_sequences(&enters, SideRole::Enters, fasta_map, flank)?;

    // Step 4. Wrap the junction homology data. The flanks are read in the call's orientation, and
    // a junction spelled from the other strand reads its GT..AG intron as CT..AC.
    let junction_homology: JunctionHomology = JunctionHomology {
        h_5prime: common_suffix(&a_in, &b_out),
        h_3prime: common_prefix(&a_out, &b_in),
        canonical_splice: (a_out.starts_with("GT") && b_out.ends_with("AG"))
            || (a_out.starts_with("CT") && b_out.ends_with("AC"))
    };

    Some(junction_homology)
}


/// Fetch the side reference sequences for the given junction side.
///
/// # Returns
/// Option<(`inside` sequence, `beyond` sequence)>
/// where `inside` sequence is the sequence the read aligned to this side junction of the junction
/// and `beyond` sequence is the template's own continuation, which the read did not follow.
fn get_side_flank_sequences(
    side: &JunctionSide,
    role: SideRole,
    fasta_map: &FastaMap,
    flank: u32
) -> Option<(String, String)> {
    let p: ReferencePosition = side.position;
    if p == 0 {
        return None;
    }

    let get = |start: ReferencePosition, end: ReferencePosition| -> Option<Box<str>> {
        if start > end {
            return Some("".into());
        }
        fasta_map
            .try_get_sequence(side.chromosome, start as usize, end as usize)
            .map(Into::into)
    };

    let (inside, beyond) = match (role, &side.strand) {
        (SideRole::Leaves, Strand::Forward) => (
            // Read moves genomically rightward and leaves: aligned tail is [p-f+1, p],
            // the continuation it abandoned is [p+1, p+f].
            get(p.saturating_sub(flank - 1).max(1), p)?.to_string(),
            get(p + 1, p + flank)?.to_string()
        ),
        (SideRole::Enters, Strand::Forward) => (
            // Read moves rightward and enters: segment head is [p, p+f-1]; what preceded it
            // on this template is [p-f, p-1].
            get(p, p + flank - 1)?.to_string(),
            get(p.saturating_sub(flank).max(1), p - 1)?.to_string()
        ),
        (SideRole::Leaves, Strand::Reverse) => (
            // Read moves genomically leftward and leaves: aligned tail is [p, p+f-1] read as
            // reverse complement; the abandoned continuation is [p-f, p-1] likewise.
            reverse_complement(&get(p, p + flank - 1)?).to_string(),
            reverse_complement(&get(p.saturating_sub(flank).max(1), p - 1)?).to_string()
        ),
        (SideRole::Enters, Strand::Reverse) => (
            // Read moves leftward and enters: segment head is [p-f+1, p] reverse complemented;
            // the preceding template sequence is [p+1, p+f] likewise.
            reverse_complement(&get(p.saturating_sub(flank - 1).max(1), p)?).to_string(),
            reverse_complement(&get(p + 1, p + flank)?).to_string()
        ),
        _ => {
            // Unreachable in practice: side_role returns None for Both/Unknown strands, so a
            // role is only ever paired with a stranded side.
            return None
        }
    };
    Some((inside.to_uppercase(), beyond.to_uppercase()))
}


fn is_foldback_shape(operation: &GraphOperation) -> bool {
    operation.get_chromosome_1() == operation.get_chromosome_2()
        && operation.get_operation_type_1() == operation.get_operation_type_2()
        && matches!(
            operation.get_operation_type_1(),
            GraphOperationType::Downstream | GraphOperationType::Upstream
        )
        && operation.get_strand_1() != operation.get_strand_2()
}


fn is_junction_variant_type(variant_type: &VariantType) -> bool {
    matches!(
        variant_type,
        VariantType::Breakpoint
            | VariantType::Translocation
            | VariantType::FusionGene
            | VariantType::CircularRNA
            | VariantType::NonCanonicalSplicing
    )
}


fn pooled_breakpoint_spread(variant_call: &VariantCall) -> Option<u32> {
    if !matches!(
        variant_call.get_consensus_graph_operation().get_variant_type(),
        VariantType::Breakpoint | VariantType::Translocation
    ) {
        return None;
    }
    // Unresolved records (Noop second side) carry no meaningful position 2; the consensus
    // builder skips them for the same reason.
    let resolved: Vec<&VariantRecord> = variant_call
        .get_variant_records()
        .iter()
        .filter(|record| *record.get_operation_2() != GraphOperationType::Noop)
        .collect();
    if resolved.len() < 2 {
        return None;
    }
    let spread = |position: fn(&VariantRecord) -> ReferencePosition| -> u32 {
        let (min, max): (ReferencePosition, ReferencePosition) = resolved
            .iter()
            .fold((u32::MAX, u32::MIN), |(min, max), record| {
                let p: ReferencePosition = position(record);
                (min.min(p), max.max(p))
            });
        max - min
    };
    Some(spread(VariantRecord::get_position_1).max(spread(VariantRecord::get_position_2)))
}


pub(crate) fn read_leaves_footprint(
    transcript_model: &TranscriptModel,
    operation: &GraphOperation,
    slack: u32
) -> bool {
    let chromosome: ReferenceChromosomeID = operation.get_chromosome_1();
    // Reference interval covered by this read's exons on one strand of the fold's
    // chromosome; `None` when there are none.
    let span_on = |strand: &Strand| -> Option<(ReferencePosition, ReferencePosition)> {
        transcript_model
            .get_exons()
            .iter()
            .filter(|exon| exon.reference_chromosome_id == chromosome
                && exon.reference_strand == *strand)
            .fold(None, |span, exon| match span {
                None => Some((exon.reference_start, exon.reference_end)),
                Some((start, end)) => Some((
                    start.min(exon.reference_start),
                    end.max(exon.reference_end)
                ))
            })
    };
    let (Some(a), Some(b)) = (
        span_on(operation.get_strand_1()),
        span_on(operation.get_strand_2())
    ) else {
        return false;   // an arm with no exon here has nothing to leave
    };
    let contains = |outer: (ReferencePosition, ReferencePosition), inner: (ReferencePosition, ReferencePosition)| -> bool {
        inner.0 + slack >= outer.0 && inner.1 <= outer.1 + slack
    };
    !(contains(a, b) || contains(b, a))
}


pub fn is_template_switch(
    verdict: TemplateSwitchVerdict,
    variant_type: &VariantType
) -> bool {
    matches!(variant_type, VariantType::Breakpoint | VariantType::Translocation)
        && matches!(verdict, TemplateSwitchVerdict::Flagged(reason) if reason != TemplateSwitchReason::SoftHomologyAtHub)
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/template_switch.rs"]
mod tests;