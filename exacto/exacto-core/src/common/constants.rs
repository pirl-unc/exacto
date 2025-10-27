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


use phf::phf_map;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::str::FromStr;


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

#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AnalyteType {
    DNA,
    RNA
}

impl AnalyteType {
    pub fn as_str(&self) -> &str {
        match self {
            AnalyteType::DNA => "dna",
            AnalyteType::RNA => "rna"
        }
    }
}

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

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum GenicRegion {
    Exonic,
    Intergenic,
    Intronic
}

impl GenicRegion {
    pub fn as_str(&self) -> &str {
        match self {
            GenicRegion::Exonic => "exonic",
            GenicRegion::Intergenic => "intergenic",
            GenicRegion::Intronic => "intronic"
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum Nucleotide {
    T,
    C,
    G,
    A,
    U,
    N,
    t,
    c,
    g,
    a,
    u,
    n
}

impl Nucleotide {
    pub fn as_str(&self) -> &str {
        match self {
            Nucleotide::T => "T",
            Nucleotide::C => "C",
            Nucleotide::G => "G",
            Nucleotide::A => "A",
            Nucleotide::U => "U",
            Nucleotide::N => "N",
            Nucleotide::t => "t",
            Nucleotide::c => "c",
            Nucleotide::g => "g",
            Nucleotide::a => "a",
            Nucleotide::u => "u",
            Nucleotide::n => "n"
        }
    }
    
    pub fn complement(&self) -> Nucleotide {
        match self {
            Nucleotide::T => Nucleotide::A,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::A => Nucleotide::T,
            Nucleotide::U => Nucleotide::A,
            Nucleotide::N => Nucleotide::N,
            Nucleotide::t => Nucleotide::a,
            Nucleotide::c => Nucleotide::g,
            Nucleotide::g => Nucleotide::c,
            Nucleotide::a => Nucleotide::t,
            Nucleotide::u => Nucleotide::a,
            Nucleotide::n => Nucleotide::n
        }
    }
}

impl FromStr for Nucleotide {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "T" => Ok(Nucleotide::T),
            "C" => Ok(Nucleotide::C),
            "G" => Ok(Nucleotide::G),
            "A" => Ok(Nucleotide::A),
            "U" => Ok(Nucleotide::U),
            "N" => Ok(Nucleotide::N),
            "t" => Ok(Nucleotide::t),
            "c" => Ok(Nucleotide::c),
            "g" => Ok(Nucleotide::g),
            "a" => Ok(Nucleotide::a),
            "u" => Ok(Nucleotide::u),
            "n" => Ok(Nucleotide::n),
            _ => Err(()),
        }
    }
}

#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
    Unknown
}

impl Strand {
    pub fn as_str(&self) -> &str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Both => "*",
            Strand::Unknown => ""
        }
    }
}

impl FromStr for Strand {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "*" => Ok(Strand::Both),
            "" => Ok(Strand::Unknown),
            _ => Err(()),
        }
    }
}
