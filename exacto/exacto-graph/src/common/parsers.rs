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


use crate::common::constants::*;
use crate::structs::vargraph_variant_node::{VarGraphVariantNode};


// pub fn parse_graph_operation(graph_operation: &str) -> VarGraphVariantNode {
//     assert!(!graph_operation.is_empty(), "Empty sequence operation.");
//     let elements: Vec<&str> = graph_operation.split(':').collect();
//     let chromosome_1: Box<str> = elements[0].into();
//     let position_1: usize = elements[1].parse().unwrap();
//     let strand_1: VarGraphStrands;
//     if elements[2] == VarGraphStrands::Forward.as_str() {
//         strand_1 = VarGraphStrands::Forward;
//     } else if elements[2] == VarGraphStrands::Reverse.as_str(){
//         strand_1 = VarGraphStrands::Reverse;
//     } else {
//         panic!("Unknown strand_1: {}", elements[2]);
//     }
//     let orientation_1: VarGraphOrientations;
//     if elements[3] == VarGraphOrientations::Upstream.as_str() {
//         orientation_1 = VarGraphOrientations::Upstream;
//     } else if elements[3] == VarGraphOrientations::Downstream.as_str(){
//         orientation_1 = VarGraphOrientations::Downstream;
//     } else {
//         panic!("Unknown orientation_1: {}", elements[3]);
//     }
//     let chromosome_2: Box<str> = elements[4].into();
//     let position_2: usize = elements[5].parse().unwrap();
//     let strand_2: VarGraphStrands;
//     if elements[6] == VarGraphStrands::Forward.as_str() {
//         strand_2 = VarGraphStrands::Forward;
//     } else if elements[6] == VarGraphStrands::Reverse.as_str(){
//         strand_2 = VarGraphStrands::Reverse;
//     } else {
//         panic!("Unknown strand_2: {}", elements[6]);
//     }
//     let orientation_2: VarGraphOrientations;
//     if elements[7] == VarGraphOrientations::Upstream.as_str() {
//         orientation_2 = VarGraphOrientations::Upstream;
//     } else if elements[7] == VarGraphOrientations::Downstream.as_str(){
//         orientation_2 = VarGraphOrientations::Downstream;
//     } else {
//         panic!("Unknown orientation_2: {}", elements[7]);
//     }
//     let sequence: Box<str> = elements[8].into();
//     let sequence_length: usize = elements[9].parse().unwrap();
//     assert!(sequence_length == sequence.len(), "Sequence lengths mismatch.");
//     let variant_node: VarGraphVariantNode = VarGraphVariantNode::new(
//         &*chromosome_1,
//         position_1,
//         strand_1,
//         orientation_1,
//         &*chromosome_2,
//         position_2,
//         strand_2,
//         orientation_2,
//         &*sequence
//     );
//     variant_node
// }
