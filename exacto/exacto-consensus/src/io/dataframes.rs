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


use polars::prelude::*;

use crate::prelude::*;


pub fn consensus_sequence_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = ConsensusSequenceRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut consensus_sequence: Vec<String> = Vec::new();
    let mut num_reads: Vec<u32> = Vec::new();
    let mut read_names: Vec<String> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        consensus_sequence.push(r.consensus_sequence.into());
        num_reads.push(r.num_reads as u32);
        read_names.push(r.read_names.into());
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("consensus_sequence".into(), consensus_sequence)),
        Column::from(Series::new("num_reads".into(), num_reads)),
        Column::from(Series::new("read_names".into(), read_names))
    ]).unwrap()
}
