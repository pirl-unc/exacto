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


extern crate bam;

use polars::prelude::*;
use polars::df;
use rustc_hash::FxHashMap;


pub fn copy_i32_series_as_vector(s: &Series) -> Vec<i32> {
    let mut v: Vec<i32> = Vec::new();
    for i in s.iter() {
        v.push(i.try_extract().unwrap());
    }
    return v;
}

pub fn copy_string_series_as_vector(s: &Series) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for i in s.iter() {
        v.push(i.to_string());
    }
    return v;
}

pub fn get_chromosome_names(bam_file: &str) -> Vec<String> {
    let reader = bam::IndexedReader::from_path(bam_file).unwrap();
    let mut chromosomes: Vec<String> = Vec::new();
    for curr_chr in reader.header().reference_names().iter() {
        chromosomes.push(curr_chr.to_string());
    }
    return chromosomes;
}