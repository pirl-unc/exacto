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
pub enum SequenceOperationType {
    Downstream,
    Read,
    Upstream,
    Mark,
    Skip
}

impl SequenceOperationType {
    pub fn as_str(&self) -> &str {
        match self {
            SequenceOperationType::Downstream => "D",
            SequenceOperationType::Read => "R",
            SequenceOperationType::Upstream => "U",
            SequenceOperationType::Mark => "M",
            SequenceOperationType::Skip => "S"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum Strand {
    Forward,
    Reverse
}

impl Strand {
    pub fn as_str(&self) -> &str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum SequenceOperationVariantType {
    Alternative3PrimeSpliceSite,
    Alternative5PrimeSpliceSite,
    Breakpoint,
    CrypticExon,
    Deletion,
    ExonSkipping,
    FusionGene,
    Insertion,
    IntronRetention,
    MultiNucleotideVariant,
    SingleNucleotideVariant,
    Translocation
}

impl SequenceOperationVariantType {
    pub fn as_str(&self) -> &str {
        match self {
            SequenceOperationVariantType::Alternative3PrimeSpliceSite => "A3P",
            SequenceOperationVariantType::Alternative5PrimeSpliceSite => "A5P",
            SequenceOperationVariantType::Breakpoint => "BND",
            SequenceOperationVariantType::CrypticExon => "CRX",
            SequenceOperationVariantType::Deletion => "DEL",
            SequenceOperationVariantType::ExonSkipping => "SKP",
            SequenceOperationVariantType::FusionGene => "FUS",
            SequenceOperationVariantType::Insertion => "A5P",
            SequenceOperationVariantType::IntronRetention => "IRT",
            SequenceOperationVariantType::MultiNucleotideVariant => "MNV",
            SequenceOperationVariantType::SingleNucleotideVariant => "SNV",
            SequenceOperationVariantType::Translocation => "TRA"
        }
    }
}
