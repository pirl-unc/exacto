use std::str::FromStr;
use serde::{Deserialize, Serialize};


/// AnalyteType
#[derive(Clone,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum AnalyteType {
    CDNA,
    DNA,
    RNA,
}

impl AnalyteType {
    pub fn as_str(&self) -> &str {
        match self {
            AnalyteType::CDNA => "cdna",
            AnalyteType::DNA => "dna",
            AnalyteType::RNA => "rna"
        }
    }
}


/// GenicRegion
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


/// Nucleotide
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
        match s {
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
            _ => Err(())
        }
    }
}


/// Strand
#[repr(u8)]
#[derive(Clone,Debug,Eq,Hash,Ord,PartialEq,PartialOrd,Serialize,Deserialize)]
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
            // GTF and BED write "." for a feature with no strand.
            "" | "." => Ok(Strand::Unknown),
            _ => Err(()),
        }
    }
}


/// TranscriptTerminus
#[repr(u8)]
#[derive(Clone,Copy,Debug,Eq,Hash,PartialEq,Serialize,Deserialize)]
pub enum TranscriptTerminus {
    FivePrime,
    ThreePrime
}


/// TranslationStrategy
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
