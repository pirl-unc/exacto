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


use exacto_caller::prelude::{GraphOperation, VariantType};
use std::collections::{HashMap, HashSet};


#[derive(Default)]
pub struct TrustedAlleles {
    /// HashMap<(chromosome ID, reference position), allele>
    pub substitutions: HashMap<(u16, u32), u8>,

    /// HashMap<(chromosome ID, reference-left insertion anchor), insertion sequence>
    pub insertions: HashMap<(u16, u32), Box<str>>,

    /// Vec<(chromosome ID, start, end)>
    pub deletions: Vec<(u16, u32, u32)>
}

impl TrustedAlleles {
    pub fn new(operations: &HashSet<&GraphOperation>) -> Self {
        let mut result: TrustedAlleles = Self::default();
        for operation in operations {
            if !matches!(
                operation.get_variant_type(),
                VariantType::SingleNucleotideVariant
                    | VariantType::MultiNucleotideVariant
                    | VariantType::Insertion
                    | VariantType::Deletion
            ) {
                // Structural operations are outside this small-variant planner.
                continue;
            }

            let chromosome: u16 = operation.get_chromosome_1();
            let left: u32 = operation.get_position_1();
            let right: u32 = operation.get_position_2();
            let sequence: String = operation.get_standardized_sequence();

            match operation.get_variant_type() {
                VariantType::SingleNucleotideVariant
                | VariantType::MultiNucleotideVariant => {
                    for (offset, allele) in sequence.bytes().enumerate() {
                        let key: (u16, u32) = (chromosome, left + 1 + offset as u32);
                        result.substitutions.insert(key, allele);
                    }
                },
                VariantType::Insertion => {
                    let key: (u16, u32) = (chromosome, left);
                    result.insertions.insert(key, sequence.into());
                },
                VariantType::Deletion => {
                    let deletion: (u16, u32, u32) = (chromosome, left, right);
                    result.deletions.push(deletion);
                }

                // The guard at the top of the loop admits only the four kinds above.
                _ => unreachable!("structural operations were skipped above")
            }
        }

        result
    }

    pub fn deletes(&self, chromosome: u16, position: u32) -> bool {
        self.deletions.iter().any(|&(chr, left, right)| {
            chr == chromosome && left < position && position < right
        })
    }
}