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


use exacto_core::prelude::*;
use std::str::FromStr;


use crate::prelude::{
    AlignmentModelBaseKind, 
    AlignmentModelRecordType, 
    AssembledTranscriptModelAlignmentRecord
};


#[derive(Clone,Debug)]
pub struct TerminalAnchor {
    pub chromosome: ReferenceChromosomeName,
    pub position: ReferencePosition,
    pub strand: Strand,
    
    /// Read position of the terminal aligned base: `read_start` of the anchor row at
    /// the 5' end, `read_end` at the 3' end.
    pub read_position: ReadPosition,

    pub reference_gene_id: ReferenceGeneID,
    pub reference_gene_name: ReferenceGeneName,
    pub reference_transcript_id: ReferenceTranscriptID
}


/// Locates the terminal reference base of an assembled transcript from its model alignment rows.
pub fn locate_terminal_anchor(
    rows: &[&AssembledTranscriptModelAlignmentRecord],
    terminus: TranscriptTerminus
) -> Option<TerminalAnchor> {
    // Step 1. Order the base rows terminal-first.
    let mut base_rows: Vec<&AssembledTranscriptModelAlignmentRecord> = rows.iter()
        .filter(|row| &*row.record_type == AlignmentModelRecordType::Base.as_str())
        .map(|row| *row)
        .collect();
    match terminus {
        TranscriptTerminus::FivePrime => base_rows.sort_by(|a, b| a.index.cmp(&b.index)),
        TranscriptTerminus::ThreePrime => base_rows.sort_by(|a, b| b.index.cmp(&a.index))
    }
    if base_rows.is_empty() {
        return None;
    }

    let strand: Strand = Strand::from_str(&base_rows[0].strand_1).unwrap_or(Strand::Unknown);
    if strand != Strand::Forward && strand != Strand::Reverse {
        return None;
    }

    // Step 2. Walk inward to the nearest row that consumed reference sequence.
    let anchor_row: &AssembledTranscriptModelAlignmentRecord = base_rows.iter()
        .find(|row| &*row.kind == AlignmentModelBaseKind::Match.as_str()
            || &*row.kind == AlignmentModelBaseKind::Mismatch.as_str())
        .map(|row| *row)?;

    // Step 3. Decode the terminal base coordinate.
    let is_mismatch: bool = &*anchor_row.kind == AlignmentModelBaseKind::Mismatch.as_str();
    let use_position_2: bool = match terminus {
        TranscriptTerminus::FivePrime => strand == Strand::Reverse,
        TranscriptTerminus::ThreePrime => strand == Strand::Forward
    };
    let (position, reference_gene_id, reference_transcript_id) = if use_position_2 {
        (if is_mismatch { anchor_row.position_2.saturating_sub(1) } else { anchor_row.position_2 },
         anchor_row.reference_gene_id_2.clone(),
         anchor_row.reference_transcript_id_2.clone())
    } else {
        (if is_mismatch { anchor_row.position_1 + 1 } else { anchor_row.position_1 },
         anchor_row.reference_gene_id_1.clone(),
         anchor_row.reference_transcript_id_1.clone())
    };
    let read_position: ReadPosition = match terminus {
        TranscriptTerminus::FivePrime => anchor_row.read_start,
        TranscriptTerminus::ThreePrime => anchor_row.read_end
    };

    Some(TerminalAnchor {
        chromosome: anchor_row.chromosome_1.clone(),
        position: position,
        strand: strand,
        read_position: read_position,
        reference_gene_id: reference_gene_id,
        reference_gene_name: resolve_reference_gene_name(anchor_row, &reference_transcript_id),
        reference_transcript_id: reference_transcript_id
    })
}


/// Converts a reference position to a read position.
pub fn reference_position_to_read_position(
    rows: &[&AssembledTranscriptModelAlignmentRecord],
    chromosome: &str,
    position: ReferencePosition,
    strand: &Strand,
    terminus: TranscriptTerminus
) -> Option<ReadPosition> {
    let mut read_position: Option<ReadPosition> = None;
    for row in rows.iter() {
        if &*row.record_type != AlignmentModelRecordType::Base.as_str()
            || &*row.kind != AlignmentModelBaseKind::Match.as_str() {
            continue;
        }
        if &*row.chromosome_1 != chromosome 
            || position < row.position_1 
            || position > row.position_2 {
            continue;
        }
        // `position_1` pairs with `read_start` on the forward strand and with
        // `read_end` on the reverse strand.
        let candidate: ReadPosition = if *strand == Strand::Forward {
            row.read_start + (position - row.position_1)
        } else {
            row.read_end - (position - row.position_1)
        };
        let keep: bool = match read_position {
            None => true,
            Some(current) => match terminus {
                TranscriptTerminus::FivePrime => candidate < current,
                TranscriptTerminus::ThreePrime => candidate > current
            }
        };
        if keep {
            read_position = Some(candidate);
        }
    }
    read_position
}


/// Picks the gene name matching a reference transcript ID.
///
/// # Notes
/// * `reference_gene_name` and `reference_transcript_id` are `;`-joined lists across
///   every reference transcript match on the model, aligned position for position.
fn resolve_reference_gene_name(
    row: &AssembledTranscriptModelAlignmentRecord,
    reference_transcript_id: &str
) -> ReferenceGeneName {
    if reference_transcript_id.is_empty() {
        return "".into();
    }
    let reference_transcript_ids: Vec<&str> = row.reference_transcript_id.split(LIST_SEPARATOR).collect();
    let reference_gene_names: Vec<&str> = row.reference_gene_name.split(LIST_SEPARATOR).collect();
    for (i, transcript_id) in reference_transcript_ids.iter().enumerate() {
        if *transcript_id == reference_transcript_id && i < reference_gene_names.len() {
            return reference_gene_names[i].into();
        }
    }
    "".into()
}
