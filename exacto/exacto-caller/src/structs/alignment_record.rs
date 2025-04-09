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


use noodles_bam as bam;


#[derive(Debug)]
pub struct AlignmentRecord {
    pub read_start: usize,
    pub read_end: usize,
    pub reverse_complemented: bool,
    pub record: bam::Record
}

impl AlignmentRecord {
    pub fn new(
        read_start: usize,
        read_end: usize,
        reverse_complemented: bool,
        record: bam::Record
    ) -> Self {
        Self {
            read_start: read_start,
            read_end: read_end,
            reverse_complemented: reverse_complemented,
            record
        }
    }
}

impl Clone for AlignmentRecord {
    fn clone(&self) -> Self {
        AlignmentRecord {
            read_start: self.read_start,
            read_end: self.read_end,
            reverse_complemented: self.reverse_complemented,
            record: self.record.clone()
        }
    }
}
