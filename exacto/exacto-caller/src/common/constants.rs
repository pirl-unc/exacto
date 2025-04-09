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


#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum SequenceOperationTypes {
    Downstream,
    Read,
    Upstream,
    Mark,
    Skip
}

impl SequenceOperationTypes {
    pub fn as_str(&self) -> &str {
        match self {
            SequenceOperationTypes::Downstream => "D",
            SequenceOperationTypes::Read => "R",
            SequenceOperationTypes::Upstream => "U",
            SequenceOperationTypes::Mark => "M",
            SequenceOperationTypes::Skip => "S"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum Strands {
    Forward,
    Reverse
}

impl Strands {
    pub fn as_str(&self) -> &str {
        match self {
            Strands::Forward => "+",
            Strands::Reverse => "-"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum SequenceOperationVariantTypes {
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

impl SequenceOperationVariantTypes {
    pub fn as_str(&self) -> &str {
        match self {
            SequenceOperationVariantTypes::Alternative3PrimeSpliceSite => "A3P",
            SequenceOperationVariantTypes::Alternative5PrimeSpliceSite => "A5P",
            SequenceOperationVariantTypes::Breakpoint => "BND",
            SequenceOperationVariantTypes::CrypticExon => "CRX",
            SequenceOperationVariantTypes::Deletion => "DEL",
            SequenceOperationVariantTypes::ExonSkipping => "SKP",
            SequenceOperationVariantTypes::FusionGene => "FUS",
            SequenceOperationVariantTypes::Insertion => "A5P",
            SequenceOperationVariantTypes::IntronRetention => "IRT",
            SequenceOperationVariantTypes::MultiNucleotideVariant => "MNV",
            SequenceOperationVariantTypes::SingleNucleotideVariant => "SNV",
            SequenceOperationVariantTypes::Translocation => "TRA"
        }
    }
}
