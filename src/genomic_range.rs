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


extern crate serde;
use serde::{Serialize, Deserialize};


#[derive(Debug, Serialize, Deserialize)]
pub struct GenomicRange {
    pub chromosome: String,
    pub start: isize,
    pub end: isize
}

impl GenomicRange {
    pub fn new(
        chromosome: String,
        start: isize,
        end: isize) -> Self {
        Self {
            chromosome: chromosome,
            start: start,
            end: end
        }
    }

    pub fn id(&self) -> String {
        let id = format!("{}:{}-{}", self.chromosome, self.start, self.end);
        return id;
    }

    pub fn overlaps(&self, chromosome: String, start: isize, end: isize) -> bool {
        // De Morgan's law on checking for non-overlapping regions
        if chromosome == self.chromosome && start <= self.end && end >= self.start {
            return true;
        } else {
            return false;
        }
    }
}

impl Clone for GenomicRange {
    fn clone(&self) -> Self {
        GenomicRange {
            chromosome: self.chromosome.clone(),
            start: self.start,
            end: self.end
        }
    }
}
