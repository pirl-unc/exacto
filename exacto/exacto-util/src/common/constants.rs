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


use std::collections::HashSet;
use once_cell::sync::Lazy;
use phf::phf_map;
use serde::{Deserialize, Serialize};

pub static CODON_TABLE: phf::Map<&'static str, &str> = phf_map! {
    "UUU" => "F", "UUC" => "F",
    "UUA" => "L", "UUG" => "L",

    "UCU" => "S", "UCC" => "S",
    "UCA" => "S", "UCG" => "S",

    "UAU" => "Y", "UAC" => "Y",
    "UAA" => "*", "UAG" => "*",

    "UGU" => "C", "UGC" => "C",
    "UGA" => "*", "UGG" => "W",

    "CUU" => "L", "CUC" => "L",
    "CUA" => "L", "CUG" => "L",

    "CCU" => "P", "CCC" => "P",
    "CCA" => "P", "CCG" => "P",

    "CAU" => "H", "CAC" => "H",
    "CAA" => "Q", "CAG" => "Q",

    "CGU" => "R", "CGC" => "R",
    "CGA" => "R", "CGG" => "R",

    "AUU" => "I", "AUC" => "I",
    "AUA" => "I", "AUG" => "M",

    "ACU" => "T", "ACC" => "T",
    "ACA" => "T", "ACG" => "T",

    "AAU" => "N", "AAC" => "N",
    "AAA" => "K", "AAG" => "K",

    "AGU" => "S", "AGC" => "S",
    "AGA" => "R", "AGG" => "R",

    "GUU" => "V", "GUC" => "V",
    "GUA" => "V", "GUG" => "V",

    "GCU" => "A", "GCC" => "A",
    "GCA" => "A", "GCG" => "A",

    "GAU" => "D", "GAC" => "D",
    "GAA" => "E", "GAG" => "E",

    "GGU" => "G", "GGC" => "G",
    "GGA" => "G", "GGG" => "G",
};

pub static START_CODONS: Lazy<HashSet<Box<str>>> = Lazy::new(|| {
    ["AUG", "GUG", "CUG", "UUG"]
        .into_iter()
        .map(|codon| codon.to_string().into_boxed_str())
        .collect()
});

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