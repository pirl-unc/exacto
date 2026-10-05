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


use csv::ReaderBuilder;

use crate::prelude::*;


pub fn load_assembled_transcript_support_records(
    tsv_file: &str
) -> Result<Vec<AssembledTranscriptSupportRecord>, TranslatorError> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .map_err(|error| TranslatorError::File { file: tsv_file.into(), reason: error.to_string().into_boxed_str() })?;
    reader
        .deserialize()
        .enumerate()
        .map(|(record, result)| result.map_err(|error| TranslatorError::Record {
            file: tsv_file.into(),
            record: record + 1,
            reason: error.to_string().into_boxed_str()
        }))
        .collect()
}
