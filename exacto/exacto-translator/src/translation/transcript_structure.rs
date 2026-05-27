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


use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptStructure {
    pub id: usize,
    pub reference_gene_names: Vec<Box<str>>,
    pub reference_transcript_ids: Vec<Box<str>>,
    pub items: Vec<TranscriptStructureItem>
}

impl TranscriptStructure {
    pub fn new(
        id: usize,
        reference_gene_names: Vec<Box<str>>,
        reference_transcript_ids: Vec<Box<str>>
    ) -> Self {
        Self {
            id: id,
            reference_gene_names: reference_gene_names,
            reference_transcript_ids: reference_transcript_ids,
            items: Vec::new()
        }
    }

    pub fn add_item(&mut self, item: TranscriptStructureItem) {
        self.items.push(item);
    }

    pub fn get_id(&self) -> usize {
        self.id
    }
}


impl Clone for TranscriptStructure {
    fn clone(&self) -> Self {
        TranscriptStructure {
            id: self.id,
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            reference_gene_names: self.reference_gene_names.clone(),
            items: self.items.clone()
        }
    }
}

impl Default for TranscriptStructure {
    fn default() -> Self {
        Self {
            id: 0,
            reference_gene_names: Vec::new(),
            reference_transcript_ids: Vec::new(),
            items: Vec::new()
        }
    }
}
