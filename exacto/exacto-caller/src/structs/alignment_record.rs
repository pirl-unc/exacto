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
use noodles_bam as bam;
use std::sync::Arc;


#[derive(Debug)]
pub struct AlignmentRecord {
    /// Read sequence start position
    pub read_start: u32,
    
    /// Read sequence end position
    pub read_end: u32,

    /// Reference strand
    pub reference_strand: Strand,

    pub record: Arc<bam::Record>
}

impl PartialEq for AlignmentRecord {
    fn eq(&self, other: &Self) -> bool {
        self.read_start == other.read_start &&
            self.read_end == other.read_end &&
            self.record == other.record
    }
}

impl AlignmentRecord {
    pub fn new(
        read_start: u32,
        read_end: u32,
        reference_strand: Strand,
        record: Arc<bam::Record>
    ) -> Self {
        assert!(read_start <= read_end);
        Self {
            read_start,
            read_end,
            reference_strand,
            record
        }
    }
}

impl Clone for AlignmentRecord {
    fn clone(&self) -> Self {
        AlignmentRecord {
            read_start: self.read_start,
            read_end: self.read_end,
            reference_strand: self.reference_strand.clone(),
            record: self.record.clone()
        }
    }
}
