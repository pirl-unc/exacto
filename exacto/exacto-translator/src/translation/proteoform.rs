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


#[derive(Clone,Debug,Serialize,Deserialize)]
pub struct Proteoform {
    pub id: u32,
    pub orf_start: u32,
    pub orf_end: u32,
    /// One letter per codon of `[orf_start, orf_end]`, the stop codon's `*` included.
    pub sequence: Box<str>
}

impl Proteoform {
    pub fn new(
        id: u32,
        orf_start: u32,
        orf_end: u32,
        sequence: Box<str>
    ) -> Self {
        assert!(orf_start < orf_end);
        assert_eq!(3 * sequence.len() as u32, orf_end - orf_start + 1);
        Self {
            id: id,
            orf_start: orf_start,
            orf_end: orf_end,
            sequence: sequence
        }
    }

    pub fn get_id(&self) -> u32 {
        self.id
    }

    pub fn get_orf_start(&self) -> u32 {
        self.orf_start
    }

    pub fn get_orf_end(&self) -> u32 {
        self.orf_end
    }

    pub fn get_length(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }
}
