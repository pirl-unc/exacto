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


use std::collections::BTreeSet;
use exacto_core::prelude::Strand;
use exacto_core::prelude::{ClusterID, ReadDepth, ReadSupport, ReferenceChromosomeID, ReferenceGeneID};
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::str::FromStr;


/// Allele
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum Allele {
    Alternate,
    Reference,
    NotCovered
}


/// AlignmentModelContext
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
#[serde(tag = "category", content = "value")]
pub enum AlignmentModelContext {
    Base(AlignmentModelBaseContext),
    Event(AlignmentModelEventContext)
}
impl AlignmentModelContext {
    pub fn as_base(&self) -> Option<&AlignmentModelBaseContext> {
        if let AlignmentModelContext::Base(ref base_ctx) = self {
            Some(base_ctx)
        } else {
            None
        }
    }

    pub fn as_event(&self) -> Option<&AlignmentModelEventContext> {
        if let AlignmentModelContext::Event(ref event_ctx) = self {
            Some(event_ctx)
        } else {
            None
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelContext::Base(base_context) => base_context.as_str(),
            AlignmentModelContext::Event(event_context) => event_context.as_str()
        }
    }
}


/// AlignmentModelBaseContext
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentModelBaseContext {
    Exonic,
    Intronic,
    Intergenic
}
impl AlignmentModelBaseContext {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelBaseContext::Exonic => "exonic",
            AlignmentModelBaseContext::Intronic => "intronic",
            AlignmentModelBaseContext::Intergenic => "intergenic"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentModelBaseContext::Exonic => ":",
            AlignmentModelBaseContext::Intronic => "$",
            AlignmentModelBaseContext::Intergenic => "",
        }
    }
}
impl FromStr for AlignmentModelBaseContext {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "exonic" | ":" => Ok(AlignmentModelBaseContext::Exonic),
            "intronic" | "$" => Ok(AlignmentModelBaseContext::Intronic),
            "intergenic" | "" => Ok(AlignmentModelBaseContext::Intergenic),
            _ => Err(())
        }
    }
}


/// AlignmentModelKind
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
#[serde(tag = "category", content = "value")]
pub enum AlignmentModelKind {
    Base(AlignmentModelBaseKind),
    Event(AlignmentModelEventKind)
}
impl AlignmentModelKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelKind::Base(base_kind) => base_kind.as_str(),
            AlignmentModelKind::Event(event_kind) => event_kind.as_str()
        }
    }
}


/// AlignmentModelBaseKind
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentModelBaseKind {
    Match,
    Mismatch,
    Insertion,
    Unaligned,
    Softclip
}
impl AlignmentModelBaseKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelBaseKind::Match => "match",
            AlignmentModelBaseKind::Mismatch => "mismatch",
            AlignmentModelBaseKind::Insertion => "insertion",
            AlignmentModelBaseKind::Unaligned => "unaligned",
            AlignmentModelBaseKind::Softclip => "softclip"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentModelBaseKind::Match => "=",
            AlignmentModelBaseKind::Mismatch => "*",
            AlignmentModelBaseKind::Insertion => "+",
            AlignmentModelBaseKind::Unaligned => "X",
            AlignmentModelBaseKind::Softclip => "S"
        }
    }
}
impl FromStr for AlignmentModelBaseKind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "match" | "=" => Ok(AlignmentModelBaseKind::Match),
            "mismatch" | "*" => Ok(AlignmentModelBaseKind::Mismatch),
            "insertion" | "+" => Ok(AlignmentModelBaseKind::Insertion),
            "unaligned" | "X" => Ok(AlignmentModelBaseKind::Unaligned),
            "softclip" | "S" => Ok(AlignmentModelBaseKind::Softclip),
            _ => Err(())
        }
    }
}


/// AlignmentModelEventContext
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentModelEventContext {
    BackSplicing,
    CanonicalSplicing,
    FusionGene,
    NonCanonicalSplicing
}
impl AlignmentModelEventContext {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelEventContext::BackSplicing => "backsplicing",
            AlignmentModelEventContext::CanonicalSplicing => "canonical",
            AlignmentModelEventContext::FusionGene => "fusion",
            AlignmentModelEventContext::NonCanonicalSplicing => "noncanonical"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentModelEventContext::BackSplicing => "/",
            AlignmentModelEventContext::CanonicalSplicing => ">",
            AlignmentModelEventContext::FusionGene => "@",
            AlignmentModelEventContext::NonCanonicalSplicing => "^"
        }
    }
}
impl FromStr for AlignmentModelEventContext {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "backsplicing" | "/" => Ok(AlignmentModelEventContext::BackSplicing),
            "canonical" | ">" => Ok(AlignmentModelEventContext::CanonicalSplicing),
            "fusion" | "@" => Ok(AlignmentModelEventContext::FusionGene),
            "noncanonical" | "^" => Ok(AlignmentModelEventContext::NonCanonicalSplicing),
            _ => Err(())
        }
    }
}


/// AlignmentModelEventKind
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentModelEventKind {
    Breakpoint,
    Deletion,
    Splicing,
    Boundary
}
impl AlignmentModelEventKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelEventKind::Breakpoint => "breakpoint",
            AlignmentModelEventKind::Deletion => "deletion",
            AlignmentModelEventKind::Splicing => "splicing",
            AlignmentModelEventKind::Boundary => "boundary"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentModelEventKind::Breakpoint => "#",
            AlignmentModelEventKind::Deletion => "-",
            AlignmentModelEventKind::Splicing => "~",
            AlignmentModelEventKind::Boundary => "|"
        }
    }
}
impl FromStr for AlignmentModelEventKind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "breakpoint" | "#" => Ok(AlignmentModelEventKind::Breakpoint),
            "deletion" | "-" => Ok(AlignmentModelEventKind::Deletion),
            "splicing" | "~" => Ok(AlignmentModelEventKind::Splicing),
            "boundary" | "|" => Ok(AlignmentModelEventKind::Boundary),
            _ => Err(())
        }
    }
}


/// AlignmentModelRecordType
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentModelRecordType {
    Base,
    Event
}
impl AlignmentModelRecordType {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentModelRecordType::Base => "base",
            AlignmentModelRecordType::Event => "event"
        }
    }
}
impl FromStr for AlignmentModelRecordType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "base" => Ok(AlignmentModelRecordType::Base),
            "event" => Ok(AlignmentModelRecordType::Event),
            _ => Err(())
        }
    }
}


/// CSTagKind
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum CSTagKind {
    Match,
    Mismatch,
    Deletion,
    Insertion,
    Splicing
}

impl CSTagKind {
    pub fn as_str(&self) -> &str {
        match self {
            CSTagKind::Match => ":",
            CSTagKind::Mismatch => "*",
            CSTagKind::Deletion => "-",
            CSTagKind::Insertion => "+",
            CSTagKind::Splicing => "~"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            CSTagKind::Match => ":",
            CSTagKind::Mismatch => "*",
            CSTagKind::Deletion => "-",
            CSTagKind::Insertion => "+",
            CSTagKind::Splicing => "~"
        }
    }

    pub fn from_symbol_str(s: &str) -> Result<Self, ()> {
        match s {
            ":" => Ok(Self::Match),
            "*" => Ok(Self::Mismatch),
            "-" => Ok(Self::Deletion),
            "+" => Ok(Self::Insertion),
            "~" => Ok(Self::Splicing),
            _ => Err(()),
        }
    }
}

impl FromStr for CSTagKind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            ":" => Ok(Self::Match),
            "*" => Ok(Self::Mismatch),
            "-" => Ok(Self::Deletion),
            "+" => Ok(Self::Insertion),
            "~" => Ok(Self::Splicing),
            _ => Err(())
        }
    }
}


/// DNAVariantOrigin
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum DNAVariantOrigin {
    Germline,
    Somatic
}
impl DNAVariantOrigin {
    pub fn as_str(&self) -> &str {
        match self {
            DNAVariantOrigin::Germline => "germline",
            DNAVariantOrigin::Somatic => "somatic"
        }
    }
}


/// GraphOperationConsensusMethod
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum GraphOperationConsensusMethod {
    Modal,
    Constructed
}
impl GraphOperationConsensusMethod {
    pub fn as_str(&self) -> &str {
        match self {
            GraphOperationConsensusMethod::Modal => "modal",
            GraphOperationConsensusMethod::Constructed => "constructed",
        }
    }
}


/// GraphOperationType
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum GraphOperationType {
    Downstream,
    Include,
    Mark,
    Skip,
    Upstream,
    Noop
}
impl GraphOperationType {
    pub fn as_str(&self) -> &str {
        match self {
            GraphOperationType::Downstream => "D",
            GraphOperationType::Include => "I",
            GraphOperationType::Mark => "M",
            GraphOperationType::Skip => "S",
            GraphOperationType::Upstream => "U",
            GraphOperationType::Noop => "N"
        }
    }

    /// Orientation of a junction departing from an anchored alignment base.
    pub fn for_breakpoint(strand: &Strand, at_segment_start: bool) -> GraphOperationType {
        match (strand, at_segment_start) {
            (Strand::Forward, true)  => GraphOperationType::Upstream,
            (Strand::Forward, false) => GraphOperationType::Downstream,
            (Strand::Reverse, true)  => GraphOperationType::Downstream,
            (Strand::Reverse, false) => GraphOperationType::Upstream,
            _ => panic!("Unexpected strand: {}", strand.as_str())
        }
    }
}
impl FromStr for GraphOperationType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "D" => Ok(GraphOperationType::Downstream),
            "I" => Ok(GraphOperationType::Include),
            "M" => Ok(GraphOperationType::Mark),
            "S" => Ok(GraphOperationType::Skip),
            "U" => Ok(GraphOperationType::Upstream),
            "N" => Ok(GraphOperationType::Noop),
            _ => Err(()),
        }
    }
}


/// JunctionBoundary
#[derive(Clone,Copy,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum JunctionBoundary {
    /// `position_1`, the first intron base in transcript order. The insertion precedes the intron.
    Donor,

    /// `position_2`, the last intron base in transcript order. The insertion follows the intron.
    Acceptor
}


#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum NonsenseMediatedDecayVerdict {
    /// The stop codon lies more than `nmd_distance_threshold` nt upstream of the last junction.
    Predicted { distance_to_last_junction: u32 },

    /// `None`: no junction downstream of the stop codon (unspliced, or stop in the last exon).
    /// Other triggers (a long 3' UTR, an upstream ORF) are not assessed.
    NotPredicted { distance_to_last_junction: Option<u32> }
}


/// TemplateSwitchReason
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub enum TemplateSwitchReason {
    /// Fold-back with an inverted repeat long enough to prime it (panel B).
    Foldback,
    /// Junction-spanning homology at or above `min_homology` (panel A).
    StrongHomology,
    /// Soft-band homology whose pooled members spread beyond the homology interval.
    SoftHomologyWithDispersion,
    /// Soft-band homology at a breakpoint partnered with many loci. Assigned at
    /// aggregation; reported, never suppressed.
    SoftHomologyAtHub
}


/// TemplateSwitchVerdict
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub enum TemplateSwitchVerdict {
    /// Not a two-sided junction, or the flanks could not be oriented.
    NotAssessed,
    Clear,
    Flagged(TemplateSwitchReason)
}


/// UnsplicedClusterKey
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum UnsplicedClusterKey {
    /// Overlaps one or more annotated gene loci on this chromosome.
    Genes(ReferenceChromosomeID, BTreeSet<ReferenceGeneID>),
    /// Overlaps no annotated gene: a single-linkage locus component on this chromosome.
    Locus(ReferenceChromosomeID, ClusterID),
}


/// VariantCallFailure
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum VariantCallFailure {
    LowTotalDepth {
        total_depth: ReadDepth,
        min_total_depth: ReadDepth
    },
    TooFewReads {
        num_reads: ReadSupport,
        min_reads: ReadSupport,
        repeat_len: Option<u32>
    },
    NotSomatic {
        num_case_reads: ReadSupport,
        num_control_reads: ReadSupport
    },
    TemplateSwitch(TemplateSwitchReason),
    /// The call matches no variant of the list the run was restricted to.
    NotListed,
    LowAltAlleleFraction,
    /// The alternate reads' strands differ from the reference reads' at a breakpoint.
    StrandBias
}
impl VariantCallFailure {
    /// The TSV spelling.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::LowTotalDepth { .. } => "low_total_depth",
            Self::TooFewReads { .. } => "too_few_reads",
            Self::TemplateSwitch(TemplateSwitchReason::Foldback) => "template_switch_foldback",
            Self::TemplateSwitch(TemplateSwitchReason::StrongHomology) => "template_switch_strong_homology",
            Self::TemplateSwitch(TemplateSwitchReason::SoftHomologyWithDispersion) => "template_switch_soft_homology_with_dispersion",
            Self::TemplateSwitch(TemplateSwitchReason::SoftHomologyAtHub) => "template_switch_soft_homology_at_hub",
            Self::NotListed => "not_in_allowed_list",
            Self::LowAltAlleleFraction => "low_alt_allele_fraction",
            Self::StrandBias => "strand_bias",
            _ => "unknown"
        }
    }
}


/// VariantType
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum VariantType {
    Breakpoint,
    CircularRNA,
    CrypticExon,
    Deletion,
    ExonTruncation,
    FusionGene,
    Insertion,
    IntronRetention,
    MultiNucleotideVariant,
    NonCanonicalSplicing,
    SingleNucleotideVariant,
    Translocation,
    UTRExtension
}
impl VariantType {
    pub fn as_str(&self) -> &str {
        match self {
            VariantType::Breakpoint => "BND",
            VariantType::CircularRNA => "CIR",
            VariantType::CrypticExon => "CRX",
            VariantType::Deletion => "DEL",
            VariantType::ExonTruncation => "SKP",
            VariantType::FusionGene => "FUS",
            VariantType::Insertion => "INS",
            VariantType::IntronRetention => "IRT",
            VariantType::MultiNucleotideVariant => "MNV",
            VariantType::NonCanonicalSplicing => "NCS",
            VariantType::SingleNucleotideVariant => "SNV",
            VariantType::Translocation => "TRA",
            VariantType::UTRExtension => "UTR"
        }
    }
    
    pub fn is_sequence_variant(&self) -> bool {
        matches!(
            self,
            VariantType::Breakpoint 
                | VariantType::CircularRNA
                | VariantType::Deletion
                | VariantType::FusionGene
                | VariantType::Insertion
                | VariantType::MultiNucleotideVariant
                | VariantType::SingleNucleotideVariant
                | VariantType::Translocation
        )
    }
}
impl FromStr for VariantType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "BND" => Ok(Self::Breakpoint),
            "CIR" => Ok(Self::CircularRNA),
            "CRX" => Ok(Self::CrypticExon),
            "DEL" => Ok(Self::Deletion),
            "SKP" => Ok(Self::ExonTruncation),
            "FUS" => Ok(Self::FusionGene),
            "INS" => Ok(Self::Insertion),
            "IRT" => Ok(Self::IntronRetention),
            "MNV" => Ok(Self::MultiNucleotideVariant),
            "NCS" => Ok(Self::NonCanonicalSplicing),
            "SNV" => Ok(Self::SingleNucleotideVariant),
            "TRA" => Ok(Self::Translocation),
            "UTR" => Ok(Self::UTRExtension),
            _ => Err(())
        }
    }
}
