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

use crate::io::records::*;


pub fn load_integrated_variant_records(tsv_file: &str) -> Vec<IntegratedVariantRecord> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(tsv_file)
        .expect("Failed to open TSV file");
    reader
        .deserialize()
        .map(|result| result.expect("Failed to deserialize IntegratedVariantRecord row"))
        .collect()
}
