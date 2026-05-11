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


use exacto_core::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashSet;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct RNA {
    pub id: Box<str>,
    pub sequence: Box<str>
}

impl RNA {
    pub fn new(
        id: Box<str>,
        sequence: Box<str>
    ) -> Self {
        Self {
            id,
            sequence
        }
    }

    pub fn translate(&self) -> Option<Translation> {
        let mut start_codons: HashSet<&str> = HashSet::new();
        start_codons.insert("AUG");

        let peptides: Vec<(Box<str>, u32, u32, u32)> = translate(
            &*self.sequence,
            start_codons
        );

        let mut translations: Translation = Translation::new(self.clone());

        for peptide in peptides {
            let peptide: Peptide = Peptide::new(
                peptide.0,
                peptide.1,
                peptide.2,
            );
            translations.add_peptide(peptide);
        }

        if translations.get_peptides_count() == 0 {
            None
        } else {
            Some(translations)
        }
    }
}

impl Clone for RNA {
    fn clone(&self) -> Self {
        RNA {
            id: self.id.clone(),
            sequence: self.sequence.clone()
        }
    }
}
