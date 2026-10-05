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


use rayon::prelude::*;
use std::collections::HashSet;

use crate::prelude::*;


pub struct AssembledTranscriptSet {
    pub transcripts: Vec<AssembledTranscript>
}


impl AssembledTranscriptSet {
    pub fn new(transcripts: Vec<AssembledTranscript>) -> Self {
        Self { transcripts }
    }

    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, AssembledTranscript> {
        self.transcripts.iter()
    }
    
    pub fn translate(
        &mut self,
        translation_strategy: TranslationStrategy,
        start_codons: &HashSet<&str>,
        num_threads: usize
    ) -> Result<(), TranslatorError> {
        // Step 1. Translate each transcript sequence.
        // Mutate each Transcript's proteoform field in place.
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .map_err(|error| TranslatorError::ThreadPool { reason: error.to_string().into_boxed_str() })?;
        thread_pool.install(|| {
            self.transcripts.par_iter_mut().for_each(|transcript| {
                transcript.translate(&translation_strategy, start_codons);
            });
        });

        // Step 2. Assign a unique ID to each proteoform across all transcripts
        let mut proteoform_id: u32 = 1;
        for transcript in self.transcripts.iter_mut() {
            for proteoform in transcript.proteoforms.iter_mut() {
                proteoform.id = proteoform_id;
                proteoform_id += 1;
            }
        }
        Ok(())
    }
}
