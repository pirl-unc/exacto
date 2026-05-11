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
use std::str::FromStr;


#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum FrameshiftState {
    InFrame,
    FrameShifted
}

impl FrameshiftState {
    pub fn as_str(&self) -> &str {
        match self {
            FrameshiftState::InFrame => "inframe",
            FrameshiftState::FrameShifted => "frameshifted"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AminoAcidChange {
    Reference,
    Mutant
}

impl AminoAcidChange {
    pub fn as_str(&self) -> &str {
        match self {
            AminoAcidChange::Reference => "reference",
            AminoAcidChange::Mutant => "mutant"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum PrimaryStructureRecordType {
    Base,
    Event
}

impl PrimaryStructureRecordType {
    pub fn as_str(&self) -> &str {
        match self {
            PrimaryStructureRecordType::Base => "base",
            PrimaryStructureRecordType::Event => "event"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum TranslationStrategy {
    AllORFs,
    LongestORF
}

impl TranslationStrategy {
    pub fn as_str(&self) -> &str {
        match self {
            TranslationStrategy::AllORFs => "all_orfs",
            TranslationStrategy::LongestORF => "longest_orf"
        }
    }
}

impl FromStr for TranslationStrategy {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "all_orfs" => Ok(TranslationStrategy::AllORFs),
            "longest_orf" => Ok(TranslationStrategy::LongestORF),
            _ => Err(())
        }
    }
}
