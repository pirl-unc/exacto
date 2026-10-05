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


use exacto_core::prelude::LIST_SEPARATOR;
use exacto_core::prelude::ReadName;
use std::collections::{HashMap, HashSet};

use crate::prelude::*;
use crate::io::formatting::sort_reference_transcript_ids;


pub struct AssembledTranscriptRNAVariantIndex<'a> {
    records: &'a [AssembledTranscriptVariantRecord],

    /// HashMap<assembled_transcript_name< Vec<records index>>
    by_assembled_transcript_name: HashMap<ReadName, Vec<usize>>,

    /// HashMap<dna_variant_id, records index>
    by_variant_id: HashMap<u32, usize>,

    /// HashMap<(assembled_transcript_name, reference transcript ID), Vec<records index>>
    by_at_and_rt: HashMap<(ReadName, String), Vec<usize>>
}

impl<'a> AssembledTranscriptRNAVariantIndex<'a> {
    pub fn new(records: &'a [AssembledTranscriptVariantRecord]) -> Self {
        let mut by_assembled_transcript_name: HashMap<ReadName, Vec<usize>> = HashMap::new();
        let mut by_variant_id: HashMap<u32, usize> = HashMap::new();
        let mut by_at_and_rt: HashMap<(ReadName, String), Vec<usize>> = HashMap::new();

        for (i, r) in records.iter().enumerate() {
            by_assembled_transcript_name
                .entry(r.assembled_transcript_name.clone())
                .or_default()
                .push(i);

            if by_variant_id.insert(r.variant_id, i).is_some() {
                panic!(
                    "Duplicate variant_id {} in records — variant_id must be unique.",
                    r.variant_id
                );
            }

            let canonical_key: String = sort_reference_transcript_ids(
                r.reference_transcript_id.split(LIST_SEPARATOR)
            );
            by_at_and_rt
                .entry((r.assembled_transcript_name.clone(), canonical_key))
                .or_default()
                .push(i);
        }

        Self {
            records,
            by_assembled_transcript_name,
            by_variant_id,
            by_at_and_rt
        }
    }

    pub fn records(&self) -> &'a [AssembledTranscriptVariantRecord] {
        self.records
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn get_by_variant_id(&self, variant_id: u32) -> Option<&'a AssembledTranscriptVariantRecord> {
        self.by_variant_id.get(&variant_id).map(|&i| &self.records[i])
    }

    pub fn get_for_assembled_transcript_name(
        &self,
        assembled_transcript_name: &str
    ) -> Vec<&'a AssembledTranscriptVariantRecord> {
        self.by_assembled_transcript_name
            .get(assembled_transcript_name)
            .map(|indices| indices.iter().map(|&i| &self.records[i]).collect())
            .unwrap_or_default()
    }

    pub fn get_for_at_and_rt(
        &self,
        assembled_transcript_name: &str,
        reference_transcript_ids: &HashSet<&str>
    ) -> Vec<&'a AssembledTranscriptVariantRecord> {
        let canonical_key: String = sort_reference_transcript_ids(reference_transcript_ids.iter().copied());
        self.by_at_and_rt
            .get(&(assembled_transcript_name.into(), canonical_key))
            .map(|indices| indices.iter().map(|&i| &self.records[i]).collect())
            .unwrap_or_default()
    }

    pub fn assembled_transcript_names(&self) -> Vec<ReadName> {
        let mut names: Vec<ReadName> = self.by_assembled_transcript_name.keys().cloned().collect();
        names.sort_unstable();
        names
    }

    pub fn grouped_by_rt_for_assembled_transcript_name(
        &self,
        assembled_transcript_name: &str
    ) -> HashMap<&'a str, Vec<&'a AssembledTranscriptVariantRecord>> {
        let mut grouped: HashMap<&'a str, Vec<&'a AssembledTranscriptVariantRecord>> = HashMap::new();
        if let Some(indices) = self.by_assembled_transcript_name.get(assembled_transcript_name) {
            for &i in indices {
                let r = &self.records[i];
                grouped
                    .entry(r.reference_transcript_id.as_ref())
                    .or_default()
                    .push(r);
            }
        }
        grouped
    }
}
