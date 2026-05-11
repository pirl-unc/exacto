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


use serde::{Deserialize,Serialize};


#[derive(Debug,Serialize,Deserialize)]
pub struct Peptide {
    pub sequence: Box<str>,
    pub orf_start: u32,
    pub orf_end: u32
}

impl Peptide {
    pub fn new(
        sequence: Box<str>,
        orf_start: u32,
        orf_end: u32
    ) -> Self {
        assert!(&*sequence != "");
        assert!(orf_start <= orf_end);
        Self {
            sequence,
            orf_start: orf_start,
            orf_end: orf_end
        }
    }

    pub fn get_length(&self) -> u32 {
        self.sequence.len() as u32
    }
}

impl Clone for Peptide {
    fn clone(&self) -> Self {
        Peptide {
            sequence: self.sequence.clone(),
            orf_start: self.orf_start,
            orf_end: self.orf_end
        }
    }
}
