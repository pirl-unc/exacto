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
            VarGraphNodeTypes::Reference => "REFERENCE",
            VarGraphNodeTypes::Variant => "VARIANT"
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
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

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
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


#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphEdgeDirections {
    FiveToThreePrime,
    ThreeToFivePrime
}

impl VarGraphEdgeDirections {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphEdgeDirections::FiveToThreePrime => "FIVE_TO_THREE_PRIME",
            VarGraphEdgeDirections::ThreeToFivePrime => "THREE_TO_FIVE_PRIME"
        }
    }
}

#[repr(u8)]
#[derive(Debug,Clone,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]
pub enum VarGraphEdgeAttributeKeys {
    Direction
}

impl VarGraphEdgeAttributeKeys {
    pub fn as_str(&self) -> &str {
        match self {
            VarGraphEdgeAttributeKeys::Direction => "DIRECTION"
        }
    }
}
