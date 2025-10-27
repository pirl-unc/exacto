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


use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};
use std::str::FromStr;


#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
#[serde(tag = "category", content = "value")]
pub enum AlignmentStructureKind {
    Base(AlignmentStructureBaseKind),
    Event(AlignmentStructureEventKind)
}

impl AlignmentStructureKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureKind::Base(base_kind) => base_kind.as_str(),
            AlignmentStructureKind::Event(event_kind) => event_kind.as_str()
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
#[serde(tag = "category", content = "value")]
pub enum AlignmentStructureContext {
    Base(AlignmentStructureBaseContext),
    Event(AlignmentStructureEventContext)
}

impl AlignmentStructureContext {
    pub fn as_base(&self) -> Option<&AlignmentStructureBaseContext> {
        if let AlignmentStructureContext::Base(ref base_ctx) = self {
            Some(base_ctx)
        } else {
            None
        }
    }

    pub fn as_event(&self) -> Option<&AlignmentStructureEventContext> {
        if let AlignmentStructureContext::Event(ref event_ctx) = self {
            Some(event_ctx)
        } else {
            None
        }
    }
    
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureContext::Base(base_context) => base_context.as_str(),
            AlignmentStructureContext::Event(event_context) => event_context.as_str()
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentStructureBaseContext {
    Exonic,
    Intronic,
    Intergenic
}

impl AlignmentStructureBaseContext {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureBaseContext::Exonic => "exonic",
            AlignmentStructureBaseContext::Intronic => "intronic",
            AlignmentStructureBaseContext::Intergenic => "intergenic"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentStructureBaseContext::Exonic => ":",
            AlignmentStructureBaseContext::Intronic => "$",
            AlignmentStructureBaseContext::Intergenic => "",
        }
    }
}

impl FromStr for AlignmentStructureBaseContext {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "exonic" | ":" => Ok(AlignmentStructureBaseContext::Exonic),
            "intronic" | "$" => Ok(AlignmentStructureBaseContext::Intronic),
            "intergenic" | "" => Ok(AlignmentStructureBaseContext::Intergenic),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentStructureBaseKind {
    Match,
    Mismatch,
    Insertion,
    Unaligned
}

impl AlignmentStructureBaseKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureBaseKind::Match => "match",
            AlignmentStructureBaseKind::Mismatch => "mismatch",
            AlignmentStructureBaseKind::Insertion => "insertion",
            AlignmentStructureBaseKind::Unaligned => "unaligned"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentStructureBaseKind::Match => "=",
            AlignmentStructureBaseKind::Mismatch => "*",
            AlignmentStructureBaseKind::Insertion => "+",
            AlignmentStructureBaseKind::Unaligned => "X"
        }
    }
}

impl FromStr for AlignmentStructureBaseKind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "match" | "=" => Ok(AlignmentStructureBaseKind::Match),
            "mismatch" | "*" => Ok(AlignmentStructureBaseKind::Mismatch),
            "insertion" | "+" => Ok(AlignmentStructureBaseKind::Insertion),
            "unaligned" | "X" => Ok(AlignmentStructureBaseKind::Unaligned),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentStructureEventContext {
    BackSplicing,
    CanonicalSplicing,
    FusionGene,
    NonCanonicalSplicing
}

impl AlignmentStructureEventContext {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureEventContext::BackSplicing => "backsplicing",
            AlignmentStructureEventContext::CanonicalSplicing => "canonical",
            AlignmentStructureEventContext::FusionGene => "fusion",
            AlignmentStructureEventContext::NonCanonicalSplicing => "noncanonical"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentStructureEventContext::BackSplicing => "/",
            AlignmentStructureEventContext::CanonicalSplicing => ">",
            AlignmentStructureEventContext::FusionGene => "@",
            AlignmentStructureEventContext::NonCanonicalSplicing => "^"
        }
    }
}

impl FromStr for AlignmentStructureEventContext {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "backsplicing" | "/" => Ok(AlignmentStructureEventContext::BackSplicing),
            "canonical" | ">" => Ok(AlignmentStructureEventContext::CanonicalSplicing),
            "fusion" | "@" => Ok(AlignmentStructureEventContext::FusionGene),
            "noncanonical" | "^" => Ok(AlignmentStructureEventContext::NonCanonicalSplicing),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentStructureEventKind {
    Breakpoint,
    Deletion,
    Splicing,
    Boundary
}

impl AlignmentStructureEventKind {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureEventKind::Breakpoint => "breakpoint",
            AlignmentStructureEventKind::Deletion => "deletion",
            AlignmentStructureEventKind::Splicing => "splicing",
            AlignmentStructureEventKind::Boundary => "boundary"
        }
    }

    pub fn as_symbol_str(&self) -> &str {
        match self {
            AlignmentStructureEventKind::Breakpoint => "#",
            AlignmentStructureEventKind::Deletion => "-",
            AlignmentStructureEventKind::Splicing => "~",
            AlignmentStructureEventKind::Boundary => "|"
        }
    }
}

impl FromStr for AlignmentStructureEventKind {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "breakpoint" | "#" => Ok(AlignmentStructureEventKind::Breakpoint),
            "deletion" | "-" => Ok(AlignmentStructureEventKind::Deletion),
            "splicing" | "~" => Ok(AlignmentStructureEventKind::Splicing),
            "boundary" | "|" => Ok(AlignmentStructureEventKind::Boundary),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AlignmentStructureRecordType {
    Base,
    Event
}

impl AlignmentStructureRecordType {
    pub fn as_str(&self) -> &str {
        match self {
            AlignmentStructureRecordType::Base => "base",
            AlignmentStructureRecordType::Event => "event"
        }
    }
}

impl FromStr for AlignmentStructureRecordType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "base" => Ok(AlignmentStructureRecordType::Base),
            "event" => Ok(AlignmentStructureRecordType::Event),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum ReferenceTranscriptSelectionStrategy {
    TopK,
    Threshold
}

impl ReferenceTranscriptSelectionStrategy {
    pub fn as_str(&self) -> &str {
        match self {
            ReferenceTranscriptSelectionStrategy::TopK => "top_k",
            ReferenceTranscriptSelectionStrategy::Threshold => "threshold"
        }
    }
}

impl FromStr for ReferenceTranscriptSelectionStrategy {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "top_k" => Ok(ReferenceTranscriptSelectionStrategy::TopK),
            "threshold" => Ok(ReferenceTranscriptSelectionStrategy::Threshold),
            _ => Err(()),
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum ReferenceTranscriptScoringMethod {
    CosineSimilarity,      // score = cosine_similarity(transcript,reference_transcript)
    L2Distance,            // score = l2_distance(transcript,reference_transcript)
    Jaccard,               // score = num_overlaps / (num_overlaps + num_nonoverlaps)
    NetOverlap,            // score = num_overlaps - num_nonoverlaps
    Nonoverlap,            // score = num_nonoverlaps
    Overlap,               // score = num_overlaps
    WeightedNetOverlap     // score = num_overlaps - 0.5 * num_nonoverlaps
}

impl ReferenceTranscriptScoringMethod {
    pub fn as_str(&self) -> &str {
        match self {
            ReferenceTranscriptScoringMethod::CosineSimilarity => "cosine_similarity",
            ReferenceTranscriptScoringMethod::Jaccard => "jaccard",
            ReferenceTranscriptScoringMethod::L2Distance => "l2",
            ReferenceTranscriptScoringMethod::NetOverlap => "net_overlap",
            ReferenceTranscriptScoringMethod::Nonoverlap => "non_overlap",
            ReferenceTranscriptScoringMethod::Overlap => "overlap",
            ReferenceTranscriptScoringMethod::WeightedNetOverlap => "weighted_net_overlap",
        }
    }
}

impl FromStr for ReferenceTranscriptScoringMethod {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "cosine_similarity" => Ok(ReferenceTranscriptScoringMethod::CosineSimilarity),
            "jaccard" => Ok(ReferenceTranscriptScoringMethod::Jaccard),
            "l2" => Ok(ReferenceTranscriptScoringMethod::L2Distance),
            "net_overlap" => Ok(ReferenceTranscriptScoringMethod::NetOverlap),
            "nonoverlap" => Ok(ReferenceTranscriptScoringMethod::Nonoverlap),
            "overlap" => Ok(ReferenceTranscriptScoringMethod::Overlap),
            "weighted_net_overlap" => Ok(ReferenceTranscriptScoringMethod::WeightedNetOverlap),
            _ => Err(()),
        }
    }
}

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
            VariantType::SingleNucleotideVariant => "SNV",
            VariantType::Translocation => "TRA",
            VariantType::UTRExtension => "UTR"
        }
    }

    // pub fn as_symbol_str(&self) -> &str {
    //     match self {
    //         VariantType::Breakpoint => "#",
    //         VariantType::CircularRNA => "@",
    //         VariantType::CrypticExon => "?",
    //         VariantType::Deletion => "-",
    //         VariantType::ExonTruncation => "!",
    //         VariantType::FusionGene => "&",
    //         VariantType::Insertion => "+",
    //         VariantType::IntronRetention => "$",
    //         VariantType::MultiNucleotideVariant => "{",
    //         VariantType::NoncanonicalSplicing => "/",
    //         VariantType::SingleNucleotideVariant => "*",
    //         VariantType::Translocation => "^"
    //     }
    // }
    //
    // pub fn from_symbol_str(s: &str) -> Result<Self, ()> {
    //     match s {
    //         "#" => Ok(Self::Breakpoint),
    //         "@" => Ok(Self::CircularRNA),
    //         "?" => Ok(Self::CrypticExon),
    //         "-" => Ok(Self::Deletion),
    //         "!" => Ok(Self::ExonTruncation),
    //         "&" => Ok(Self::FusionGene),
    //         "+" => Ok(Self::Insertion),
    //         "$" => Ok(Self::IntronRetention),
    //         "{" => Ok(Self::MultiNucleotideVariant),
    //         "/" => Ok(Self::NoncanonicalSplicing),
    //         "*" => Ok(Self::SingleNucleotideVariant),
    //         "^" => Ok(Self::Translocation),
    //         _ => Err(())
    //     }
    // }
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
            "SNV" => Ok(Self::SingleNucleotideVariant),
            "TRA" => Ok(Self::Translocation),
            "UTR" => Ok(Self::UTRExtension),
            _ => Err(())
        }
    }
}

