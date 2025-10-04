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
use serde::{Deserialize, Serialize};
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModelIntron {
    pub reference_chromosome_id: u16,
    pub reference_start: u32,
    pub reference_end: u32,
    pub reference_strand: Strand,
    pub intron_number: u16,
    pub donor_splice_site_signal: Box<str>,
    pub acceptor_splice_site_signal: Box<str>,
    
    /// Read position immediately before the intron starts (towards 5')
    pub exon_5p_flank_read_position: u32,
    
    /// Read position immediately after the intron ends (towards 3')
    pub exon_3p_flank_read_position: u32
}

impl PartialEq for TranscriptModelIntron {
    fn eq(&self, other: &Self) -> bool {
        self.reference_chromosome_id == other.reference_chromosome_id &&
            self.reference_start == other.reference_start &&
            self.reference_end == other.reference_end &&
            self.intron_number == other.intron_number &&
            self.donor_splice_site_signal == other.donor_splice_site_signal &&
            self.acceptor_splice_site_signal == other.acceptor_splice_site_signal &&
            self.reference_strand == other.reference_strand &&
            self.exon_5p_flank_read_position == other.exon_5p_flank_read_position &&
            self.exon_3p_flank_read_position == other.exon_3p_flank_read_position
    }
}

impl Eq for TranscriptModelIntron {}

impl Hash for TranscriptModelIntron {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_chromosome_id.hash(state);
        self.reference_start.hash(state);
        self.reference_end.hash(state);
        self.intron_number.hash(state);
        self.donor_splice_site_signal.hash(state);
        self.acceptor_splice_site_signal.hash(state);
        self.reference_strand.hash(state);
        self.exon_5p_flank_read_position.hash(state);
        self.exon_3p_flank_read_position.hash(state);
    }
}

impl TranscriptModelIntron {
    pub fn new(
        reference_chromosome_id: u16,
        reference_start: u32,
        reference_end: u32,
        reference_strand: Strand,
        intron_number: u16,
        donor_splice_site_signal: &str,
        acceptor_splice_site_signal: &str,
        exon_5p_flank_read_position: u32,
        exon_3p_flank_read_position: u32
    ) -> Self {
        assert_eq!(
            exon_5p_flank_read_position, exon_3p_flank_read_position - 1,
            "exon_3p_flank_read_position should be the read position immediately following exon_5p_flank_read_position: {}-{}",
            exon_5p_flank_read_position,
            exon_3p_flank_read_position
        );
        Self {
            reference_chromosome_id,
            reference_start,
            reference_end,
            reference_strand: reference_strand.clone(),
            intron_number,
            donor_splice_site_signal: donor_splice_site_signal.into(),
            acceptor_splice_site_signal: acceptor_splice_site_signal.into(),
            exon_5p_flank_read_position,
            exon_3p_flank_read_position
        }
    }
}

impl Clone for TranscriptModelIntron {
    fn clone(&self) -> Self {
        TranscriptModelIntron {
            reference_chromosome_id: self.reference_chromosome_id,
            reference_start: self.reference_start,
            reference_end: self.reference_end,
            intron_number: self.intron_number,
            donor_splice_site_signal: self.donor_splice_site_signal.clone(),
            acceptor_splice_site_signal: self.acceptor_splice_site_signal.clone(),
            reference_strand: self.reference_strand.clone(),
            exon_5p_flank_read_position: self.exon_5p_flank_read_position,
            exon_3p_flank_read_position: self.exon_3p_flank_read_position
        }
    }
}
