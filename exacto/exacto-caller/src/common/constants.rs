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
            AlignmentStructureEventContext::CanonicalSplicing => "canonical_splicing",
            AlignmentStructureEventContext::FusionGene => "fusion_gene",
            AlignmentStructureEventContext::NonCanonicalSplicing => "non_canonical_splicing"
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
    Upstream
}

impl GraphOperationType {
    pub fn as_str(&self) -> &str {
        match self {
            GraphOperationType::Downstream => "D",
            GraphOperationType::Include => "I",
            GraphOperationType::Mark => "M",
            GraphOperationType::Skip => "S",
            GraphOperationType::Upstream => "U"
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

