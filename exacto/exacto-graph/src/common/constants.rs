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
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum MultiDiGraphDirections {
    Incoming,
    Outgoing
}

impl MultiDiGraphDirections {
    pub fn as_str(&self) -> &str {
        match self {
            MultiDiGraphDirections::Incoming => "INCOMING",
            MultiDiGraphDirections::Outgoing => "OUTGOING"
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphNodeTypes {
    Reference,
    Variant
}

impl VarGraphNodeTypes {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphNodeTypes::Reference => "REF",
            VarGraphNodeTypes::Variant => "VAR"
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphOrientations {
    Upstream,
    Downstream,
    Any
}

impl VarGraphOrientations {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphOrientations::Upstream => "U",
            VarGraphOrientations::Downstream => "D",
            VarGraphOrientations::Any => "*"
        }
    }
}

impl FromStr for VarGraphOrientations {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "U" => Ok(VarGraphOrientations::Upstream),
            "D" => Ok(VarGraphOrientations::Downstream),
            "*" => Ok(VarGraphOrientations::Any),
            _ => Err(()),
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphReferenceNodeAnnotationKeys {
    Type
}

impl VarGraphReferenceNodeAnnotationKeys {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphReferenceNodeAnnotationKeys::Type => "type"
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphReferenceNodeTypes {
    Exonic,
    Intronic,
    Intergenic
}

impl VarGraphReferenceNodeTypes {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphReferenceNodeTypes::Exonic => "exonic",
            VarGraphReferenceNodeTypes::Intronic => "intronic",
            VarGraphReferenceNodeTypes::Intergenic => "intergenic"
        }
    }
}

impl FromStr for VarGraphReferenceNodeTypes {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "exonic" => Ok(VarGraphReferenceNodeTypes::Exonic),
            "intronic" => Ok(VarGraphReferenceNodeTypes::Intronic),
            "intergenic" => Ok(VarGraphReferenceNodeTypes::Intergenic),
            _ => Err(()),
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphStrands {
    Forward,
    Reverse
}

impl VarGraphStrands {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphStrands::Forward => "+",
            VarGraphStrands::Reverse => "-"
        }
    }
}

impl FromStr for VarGraphStrands {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(VarGraphStrands::Forward),
            "-" => Ok(VarGraphStrands::Reverse),
            _ => Err(()),
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphTypes {
    Individual,
    Population
}

impl VarGraphTypes {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphTypes::Individual => "individual",
            VarGraphTypes::Population => "population"
        }
    }
}

impl FromStr for VarGraphTypes {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "individual" => Ok(VarGraphTypes::Individual),
            "population" => Ok(VarGraphTypes::Population),
            _ => Err(())
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,Hash,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphEdgeEvent {
    Enter(usize),
    Back(usize),
    Exit(usize)
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphEdgeAttributeKeys {
    Direction,
    Orientation,
    Strand
}

impl VarGraphEdgeAttributeKeys {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphEdgeAttributeKeys::Direction => "DIRECTION",
            VarGraphEdgeAttributeKeys::Orientation => "ORIENTATION",
            VarGraphEdgeAttributeKeys::Strand => "STRAND"
        }
    }
}
