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


extern crate noodles;
extern crate noodles_core;
use noodles::bam as bam;
use noodles_core::{Region, Position};
use std::collections::HashMap;


pub fn get_chromosomes(bam_file: &str) -> HashMap<usize, (String, usize)> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header();
    let mut chromosomes: HashMap<usize, (String, usize)> = HashMap::new();
    let mut i: usize = 0;
    for chromosome in header.unwrap().reference_sequences().iter() {
        let chromosome_name: &str = &chromosome.0.to_string();
        let chromosome_length: usize = chromosome.1.length().into();
        chromosomes.insert(i, (chromosome_name.to_string(), chromosome_length));
        i += 1;
    }
    return chromosomes;
}